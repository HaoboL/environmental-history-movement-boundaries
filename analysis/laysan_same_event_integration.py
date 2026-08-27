#!/usr/bin/env python3
"""D127: same-event integration of absolute CHL and last-passage placement.

No RD is rerun. Goto joins frozen D104/D43 events. Laysan reconstructs strict
radial records from frozen start-to-endpoint paths and frozen CHL cell runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
import last_record_decomposition as last_record  # noqa: E402

PREREG = ROOT / "metadata/laysan_same_event_specification_cn.md"
D104_A = ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz"
INPUTS = ROOT / "external_inputs/laysan_same_event"
D43 = INPUTS / "absolute_chl_step_event_table.csv.gz"
D41_LAYSAN = INPUTS / "laysan_standardized_events.csv.gz"
LAYSAN_PATHS = INPUTS / "laysan_paths_start_to_endpoint.csv.gz"
LAYSAN_TOKENS = INPUTS / "laysan_event_chl_cell_runs.csv.gz"
DEFAULT_OUT = ROOT / "results/laysan_same_event"
SCALES = (500, 1000, 2000)
SEED = 1270826
EARTH_R_KM = 6371.0088


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(value: Any) -> Any:
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite(v) for v in value]
    return value


class Progress:
    def __init__(self, out: Path):
        self.out = out
        self.path = out / "progress.jsonl"
        self.t0 = time.monotonic()

    def emit(self, stage: str, done: int, total: int, **extra: Any) -> None:
        elapsed = max(time.monotonic() - self.t0, 1e-9)
        rate = done / elapsed if done else 0.0
        row = {
            "time": pd.Timestamp.now(tz="America/Los_Angeles").isoformat(),
            "stage": stage, "completed": int(done), "total": int(total),
            "percent": 100 * done / max(total, 1), "elapsed_s": elapsed,
            "throughput_per_s": rate,
            "eta_s": (total - done) / rate if rate > 0 and done < total else 0.0,
            **extra,
        }
        with self.path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(finite(row), ensure_ascii=False) + "\n")
        print(f"PROGRESS {stage} {done}/{total} ({row['percent']:.1f}%) elapsed={elapsed:.1f}s eta={row['eta_s']:.1f}s", flush=True)


def unit_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    phi = np.radians(np.asarray(lat, float)); lam = np.radians(np.asarray(lon, float))
    c = np.cos(phi)
    return EARTH_R_KM * np.column_stack((c * np.cos(lam), c * np.sin(lam), np.sin(phi)))


def dense_original_indices(path: pd.DataFrame) -> np.ndarray:
    """Indices of original fixes in the exact previously used <=4-km densification."""
    xyz = unit_xyz(path.lat.to_numpy(float), path.lon.to_numpy(float)) / EARTH_R_KM
    dots = np.clip(np.sum(xyz[:-1] * xyz[1:], axis=1), -1.0, 1.0)
    distances = EARTH_R_KM * np.arccos(dots)
    intervals = np.maximum(1, np.ceil(distances / 4.0).astype(int))
    return np.r_[0, np.cumsum(intervals)].astype(int)


def reconstruct_laysan(target_ids: set[str], smoke: bool, progress: Progress) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    paths = pd.read_csv(LAYSAN_PATHS, low_memory=False)
    paths = paths[paths.delta_m.isin(SCALES) & paths.event_id.astype(str).isin(target_ids)].copy()
    tokens = pd.read_csv(LAYSAN_TOKENS, low_memory=False)
    tokens = tokens[tokens.date_shift_days.eq(0) & tokens.event_id.astype(str).isin(target_ids)].copy()
    if smoke:
        keep = set(paths.event_id.drop_duplicates().head(100).astype(str))
        paths = paths[paths.event_id.astype(str).isin(keep)]
        tokens = tokens[tokens.event_id.astype(str).isin(keep)]
    token_groups = {str(k): g.sort_values("cell_run_index") for k, g in tokens.groupby("event_id", sort=False)}
    records: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    groups = list(paths.groupby("event_id", sort=False))
    for done, (event_id, raw) in enumerate(groups, 1):
        g = raw.sort_values("path_point_index", kind="stable").reset_index(drop=True)
        tok = token_groups.get(str(event_id))
        status = "ok"
        geometry_ok = False
        n_record_runs = 0
        n_runs = 0
        if tok is None:
            status = "missing_shift0_tokens"
        else:
            values = tok.cell_run_logchl.to_numpy(float)
            counts = tok.dense_points_in_run.to_numpy(int)
            n_runs = len(values)
            dense_idx = dense_original_indices(g)
            if counts.sum() != dense_idx[-1] + 1:
                status = "densification_count_mismatch"
            else:
                xyz = unit_xyz(g.lat.to_numpy(float), g.lon.to_numpy(float))
                origin = xyz[0]
                exact_rec = []
                maximum = 0.0
                for i in range(1, len(xyz)):
                    radius = float(np.linalg.norm(xyz[i] - origin))
                    if radius > maximum:
                        maximum = radius
                        exact_rec.append(i)
                rec = np.asarray(exact_rec, int)
                geometry_ok = bool(len(rec) and rec[-1] == len(g) - 1)
                run_ends = np.cumsum(counts)
                record_runs = np.unique(np.searchsorted(run_ends, dense_idx[rec], side="right")) if len(rec) else np.asarray([], int)
                n_record_runs = len(record_runs)
                if not geometry_ok:
                    status = "endpoint_not_last_strict_record"
                elif not last_record.eligible_values(values):
                    status = "chl_not_eligible"
                elif n_record_runs < 2:
                    status = "fewer_than_2_record_runs"
                else:
                    cluster = str(g.individual.iloc[0])
                    # D41/D43 block is frozen and joined later; temporary segment is run.
                    records.append({
                        "dataset": "usgs_laysan_albatross", "scale_m": int(g.delta_m.iloc[0]),
                        "family": "A", "event_id": str(event_id), "segment_key": str(g.run_key.iloc[0]),
                        "cluster": cluster, "values": values, "position": record_runs,
                        "observed": last_record.metrics_a(values, record_runs), "rho_tau_distinct_cell": None,
                    })
        audits.append({
            "event_id": str(event_id), "scale_m": int(g.delta_m.iloc[0]), "cluster": str(g.individual.iloc[0]),
            "run_key": str(g.run_key.iloc[0]), "status": status,
            "geometry_endpoint_is_last_record": geometry_ok, "n_runs": n_runs,
            "n_record_runs": n_record_runs, "eligible": status == "ok",
        })
        if done % max(1, len(groups) // 20) == 0 or done == len(groups):
            progress.emit("reconstruct_laysan", done, len(groups), eligible=len(records))
    return records, pd.DataFrame(audits)


def within_cluster_z(frame: pd.DataFrame, column: str) -> pd.Series:
    def one(s: pd.Series) -> pd.Series:
        sd = float(s.std(ddof=0))
        return (s - s.mean()) / sd if sd > 0 else pd.Series(np.zeros(len(s)), index=s.index)
    return frame.groupby("cluster", group_keys=False)[column].transform(one)


def fit_fe(frame: pd.DataFrame) -> np.ndarray:
    y = frame["z_log_length"].to_numpy(float)
    x = frame[["z_absolute", "L_low_dm", "L_high_dm"]].to_numpy(float)
    return np.linalg.lstsq(x, y, rcond=None)[0]


def prepare_model(frame: pd.DataFrame) -> pd.DataFrame:
    d = frame.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["step_length_km", "interior_median_logchl", "L_low", "L_high", "L_union", "cluster"]
    ).copy()
    d = d[d.step_length_km > 0].copy()
    d["log_length"] = np.log(d.step_length_km)
    d["z_log_length"] = within_cluster_z(d, "log_length")
    d["z_absolute"] = within_cluster_z(d, "interior_median_logchl")
    d["L_low_dm"] = d.L_low - d.groupby("cluster").L_low.transform("mean")
    d["L_high_dm"] = d.L_high - d.groupby("cluster").L_high.transform("mean")
    for source, target in (("interior_median_logchl", "abs_tertile"), ("log_length", "length_tertile")):
        rank = d.groupby("cluster")[source].rank(method="average", pct=True)
        d[target] = np.minimum((rank * 3).astype(int), 2)
    return d


def bootstrap_summary(frame: pd.DataFrame, reps: int, rng: np.random.Generator) -> tuple[dict[str, Any], pd.DataFrame]:
    d = prepare_model(frame)
    clusters = d.cluster.drop_duplicates().to_numpy(object)
    by = {c: d[d.cluster == c] for c in clusters}
    beta = fit_fe(d)
    # Equal-weight occupied grids within animal, then equal-weight animals.
    animal_grid = d.groupby(["cluster", "abs_tertile", "length_tertile"], as_index=False).L_union.mean()
    animal_stat = animal_grid.groupby("cluster").L_union.mean()
    estimate_l = float(animal_stat.mean())
    grid = animal_grid.groupby(["abs_tertile", "length_tertile"], as_index=False).L_union.mean()
    # Cluster bootstrap through per-animal sufficient statistics. This is
    # algebraically identical to concatenating every sampled animal table, but
    # avoids 20,000 expensive pandas concatenations per model.
    xtx = np.empty((len(clusters), 3, 3), float)
    xty = np.empty((len(clusters), 3), float)
    l_values = np.empty(len(clusters), float)
    for j, c in enumerate(clusters):
        part = by[c]
        x = part[["z_absolute", "L_low_dm", "L_high_dm"]].to_numpy(float)
        y = part["z_log_length"].to_numpy(float)
        xtx[j] = x.T @ x
        xty[j] = x.T @ y
        l_values[j] = float(animal_stat.loc[c])
    choices = rng.integers(0, len(clusters), size=(reps, len(clusters)))
    counts = np.zeros((reps, len(clusters)), dtype=np.int16)
    np.add.at(counts, (np.repeat(np.arange(reps), len(clusters)), choices.ravel()), 1)
    boot_xtx = np.einsum("rc,cij->rij", counts, xtx, optimize=True)
    boot_xty = np.einsum("rc,cj->rj", counts, xty, optimize=True)
    beta_boot = np.einsum("rij,rj->ri", np.linalg.pinv(boot_xtx), boot_xty, optimize=True)
    l_boot = counts @ l_values / len(clusters)
    blo, bhi = np.quantile(beta_boot, [0.025, 0.975], axis=0)
    llo, lhi = np.quantile(l_boot, [0.025, 0.975])
    out = {
        "events": len(d), "animals": len(clusters),
        "beta_absolute": beta[0], "beta_absolute_ci_low": blo[0], "beta_absolute_ci_high": bhi[0],
        "beta_L_low": beta[1], "beta_L_low_ci_low": blo[1], "beta_L_low_ci_high": bhi[1],
        "beta_L_high": beta[2], "beta_L_high_ci_low": blo[2], "beta_L_high_ci_high": bhi[2],
        "standardized_L_union": estimate_l, "standardized_L_union_ci_low": llo, "standardized_L_union_ci_high": lhi,
        "positive_grids": int((grid.L_union > 0).sum()), "estimable_grids": len(grid), "bootstrap_reps": reps,
    }
    return out, grid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--bootstrap-reps", type=int, default=20000)
    ap.add_argument("--phase-reps", type=int, default=999)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output / "progress.jsonl").exists():
        (args.output / "progress.jsonl").unlink()
    progress = Progress(args.output)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 forbidden: {affinity}")
    if args.smoke:
        args.bootstrap_reps = min(args.bootstrap_reps, 200)
        args.phase_reps = min(args.phase_reps, 19)
    rng = np.random.default_rng(SEED)
    progress.emit("start", 0, 1, affinity=affinity, smoke=args.smoke)

    use43 = ["dataset", "scale_m", "event_id", "cluster_id", "step_length_km",
             "interior_median_logchl", "high20_excess", "low20_excess", "extreme20_excess"]
    d43 = pd.read_csv(D43, usecols=use43, low_memory=False)
    d104a = pd.read_csv(D104_A, low_memory=False)
    goto = d104a[d104a.dataset.eq("goto")].merge(
        d43[d43.dataset.eq("goto")], on=["dataset", "scale_m", "event_id"], validate="one_to_one"
    )
    goto["cluster"] = goto["cluster"].astype(str)
    if args.smoke:
        goto = goto[goto.cluster.isin(goto.cluster.drop_duplicates().head(3))].copy()
    progress.emit("join_goto", len(goto), len(goto), animals=goto.cluster.nunique())

    target = d43[d43.dataset.eq("usgs_laysan_albatross") & d43.scale_m.isin(SCALES)].copy()
    d41_keys = pd.read_csv(D41_LAYSAN, usecols=["event_id", "scale_m", "block_id"], low_memory=False)
    target = target.merge(d41_keys, on=["event_id", "scale_m"], validate="one_to_one")
    laysan_records, laysan_audit = reconstruct_laysan(set(target.event_id.astype(str)), args.smoke, progress)
    laysan_events = last_record.event_frame(laysan_records, last_record.A_NAMES)
    laysan = laysan_events.merge(target, on=["dataset", "scale_m", "event_id"], validate="one_to_one")
    laysan["cluster"] = laysan["cluster"].astype(str)
    # Use frozen D41/D43 block for the phase-null grouping.
    block_map = target.set_index(["scale_m", "event_id"])["block_id"].astype(str).to_dict()
    for rec in laysan_records:
        rec["segment_key"] = block_map[(int(rec["scale_m"]), str(rec["event_id"]))]

    # Exact E anchor against the frozen event endpoint excess.
    anchor_diffs = np.r_[
        np.abs(laysan.E_high - laysan.high20_excess),
        np.abs(laysan.E_low - laysan.low20_excess),
        np.abs(laysan.E_union - laysan.extreme20_excess),
    ]
    anchor_max = float(np.nanmax(anchor_diffs))
    identity_max = float(np.max(np.abs(laysan.E_union - laysan.R_union - laysan.L_union)))
    if anchor_max > 1e-12 or identity_max > 1e-12:
        raise RuntimeError(f"Laysan anchor/identity failed anchor={anchor_max} identity={identity_max}")
    progress.emit("anchor_laysan", 1, 1, anchor_max=anchor_max, identity_max=identity_max)

    phase_summary, phase_audit = last_record.summarize_family(
        laysan_records, "A", args.phase_reps, args.bootstrap_reps, rng, progress
    )
    rows = []
    grids = []
    for dataset, scale, frame in [("goto", 100, goto)] + [
        ("usgs_laysan_albatross", s, laysan[laysan.scale_m.eq(s)]) for s in SCALES
    ]:
        if len(frame) == 0:
            continue
        summary, grid = bootstrap_summary(frame, args.bootstrap_reps, rng)
        summary.update({"dataset": dataset, "scale_m": scale})
        rows.append(summary)
        grid.insert(0, "scale_m", scale); grid.insert(0, "dataset", dataset)
        grids.append(grid)
        progress.emit("joint_model", len(rows), 4, dataset=dataset, scale_m=scale)
    model = pd.DataFrame(rows)
    grid_table = pd.concat(grids, ignore_index=True)

    phase_l = phase_summary[phase_summary.estimand.eq("L_union")].copy()
    phase_pass = {(r.dataset, int(r.scale_m)): bool(r.positive_support) for r in phase_l.itertuples(index=False)}
    model["conditional_L_pass"] = model.standardized_L_union_ci_low.gt(0) & model.positive_grids.ge(6)
    model["absolute_pass"] = model.beta_absolute_ci_high.lt(0)
    model["L_low_positive"] = model.beta_L_low.gt(0)
    model["L_low_strict"] = model.beta_L_low_ci_low.gt(0)
    model["phase_L_union_pass"] = [phase_pass.get(("usgs_laysan_albatross", int(s)), True) if d == "usgs_laysan_albatross" else True for d, s in model[["dataset", "scale_m"]].itertuples(index=False, name=None)]
    lp = model[model.dataset.eq("usgs_laysan_albatross") & model.conditional_L_pass & model.absolute_pass & model.phase_L_union_pass]
    scales_pass = sorted(lp.scale_m.astype(int).tolist())
    adjacent = any(a in scales_pass and b in scales_pass for a, b in zip(SCALES[:-1], SCALES[1:]))
    g = model[model.dataset.eq("goto")].iloc[0]
    l_low_strict_cross = bool(g.L_low_strict and model[model.dataset.eq("usgs_laysan_albatross")].L_low_strict.any())
    verdict_pass = bool(g.conditional_L_pass and g.absolute_pass and g.L_low_positive and adjacent and l_low_strict_cross)

    model.to_csv(args.output / "joint_model_summary.csv", index=False)
    grid_table.to_csv(args.output / "conditional_L_3x3_grid.csv", index=False)
    laysan_events.to_csv(args.output / "laysan_last_passage_event_metrics.csv.gz", index=False, compression="gzip")
    laysan_audit.to_csv(args.output / "laysan_reconstruction_audit.csv.gz", index=False, compression="gzip")
    phase_summary.to_csv(args.output / "laysan_phase_summary.csv", index=False)
    phase_audit.to_csv(args.output / "laysan_phase_null_audit.csv.gz", index=False, compression="gzip")
    payload = {
        "status": "SMOKE_COMPLETE" if args.smoke else "FORMAL_COMPLETE",
        "verdict": "SAME_EVENT_DUAL_REFERENCE_NONREDUNDANCY_SUPPORTED" if verdict_pass else "SAME_EVENT_DUAL_REFERENCE_GATE_FAILED",
        "gate_pass": verdict_pass, "laysan_scales_passing_all": scales_pass,
        "laysan_adjacent_scale_pass": adjacent, "cross_system_L_low_strict": l_low_strict_cross,
        "affinity": affinity, "cpu1_excluded": 1 not in affinity,
        "anchor_max_abs_E_difference": anchor_max, "identity_max_abs_error": identity_max,
        "model_summary": model.to_dict("records"),
        "phase_L_union": phase_l.to_dict("records"),
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (D104_A, D43, D41_LAYSAN, LAYSAN_PATHS, LAYSAN_TOKENS, PREREG)},
        "script_sha256": sha256(Path(__file__)),
        "interpretation_limit": "Same-event retrospective statistical nonredundancy only; not online cue, sensory mechanism, causal control, prey success, or Levy generation.",
    }
    (args.output / "final_summary.json").write_text(json.dumps(finite(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("complete", 1, 1, verdict=payload["verdict"])


if __name__ == "__main__":
    main()
