#!/usr/bin/env python3
"""D128: macro phenotype counterfactual under alternative movement boundaries.

The frozen Goto RD100 event sequence is retained. Boundaries are the frozen
last maximum (rho), frozen confirmation fix (tau), or one uniformly sampled
strict radial record per event. No RD detection is rerun.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
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
from scipy.stats import spearmanr

ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
from paper2_core import fit_step_distributions, xy_from_latlon  # noqa: E402

PREREG = ROOT / "metadata/boundary_counterfactual_specification_cn.md"
INPUTS = ROOT / "external_inputs/boundary_counterfactual"
CANON = INPUTS / "goto_canonical_events.csv"
DENSE = INPUTS / "goto_continuous_paths_with_chl.csv.gz"
OLD = INPUTS / "observed_segment_tail_metrics.csv"
DEFAULT_OUT = ROOT / "results/boundary_counterfactual"
SEED = 1280826
_SEGMENTS: list[dict[str, Any]] = []


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(value: Any) -> Any:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite(v) for v in value]
    return value


class Progress:
    def __init__(self, out: Path):
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


def strict_records(x: np.ndarray, y: np.ndarray, start: int, rho: int) -> np.ndarray:
    maximum = 0.0; records = []
    for i in range(start + 1, rho + 1):
        r = float(math.hypot(x[i] - x[start], y[i] - y[start]))
        if r > maximum:
            maximum = r; records.append(i)
    return np.asarray(records, dtype=int)


def lengths_from_boundaries(x: np.ndarray, y: np.ndarray, initial: int, boundaries: np.ndarray) -> np.ndarray:
    idx = np.r_[int(initial), np.asarray(boundaries, int)]
    return np.hypot(np.diff(x[idx]), np.diff(y[idx])) / 1000.0


def self_tests() -> None:
    x = np.array([0., 2., 1., 3.]); y = np.zeros(4)
    if not np.array_equal(strict_records(x, y, 0, 3), np.array([1, 3])):
        raise AssertionError("strict-record self-test failed")
    lengths = lengths_from_boundaries(x * 1000, y, 0, np.array([1, 3]))
    if not np.allclose(lengths, [2., 1.]):
        raise AssertionError("boundary-length self-test failed")


def build_segments(smoke: bool, progress: Progress) -> tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(OLD, low_memory=False)
    if smoke:
        old = old.head(20).copy()
    wanted = set(old.segment_key.astype(str))
    cols = ["track_id", "segment_id", "step_id", "is_terminal_step", "start_orig_idx",
            "endpoint_orig_idx", "trigger_orig_idx", "step_length_km"]
    events = pd.read_csv(CANON, usecols=cols, low_memory=False)
    events["segment_key"] = events.track_id.astype(str) + "|" + events.segment_id.astype(str)
    events = events[events.segment_key.isin(wanted) & ~events.is_terminal_step.astype(bool)].copy()
    dense = pd.read_csv(DENSE, low_memory=False)
    dense = dense[dense.orig_idx.notna()].copy()
    dense["segment_key"] = dense.track_id.astype(str) + "|" + dense.segment_id.astype(str)
    dense = dense[dense.segment_key.isin(wanted)].copy()
    groups = {str(k): g.sort_values("orig_idx", kind="stable") for k, g in dense.groupby("segment_key", sort=False)}
    rows: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    for done, oldrow in enumerate(old.itertuples(index=False), 1):
        key = str(oldrow.segment_key)
        g = groups[key]
        orig = g.orig_idx.to_numpy(int)
        if not np.array_equal(orig, np.arange(len(g))):
            raise RuntimeError(f"noncontiguous original index: {key}")
        x, y = xy_from_latlon(g.Lat.to_numpy(float), g.Lon.to_numpy(float))
        ev = events[events.segment_key.eq(key)].sort_values("step_id", kind="stable")
        starts = ev.start_orig_idx.to_numpy(int)
        rhos = ev.endpoint_orig_idx.to_numpy(int)
        taus = ev.trigger_orig_idx.to_numpy(float)
        if not np.isfinite(taus).all():
            raise RuntimeError(f"nonterminal event lacks tau: {key}")
        taus = taus.astype(int)
        candidates = [strict_records(x, y, int(a), int(b)) for a, b in zip(starts, rhos)]
        if any(len(c) == 0 or c[-1] != r for c, r in zip(candidates, rhos)):
            raise RuntimeError(f"rho is not final strict record: {key}")
        rho_coord = np.hypot(x[rhos] - x[starts], y[rhos] - y[starts]) / 1000.0
        frozen = ev.step_length_km.to_numpy(float)
        rho_error = float(np.max(np.abs(rho_coord - frozen)))
        tau_nonincreasing = int(np.sum(np.diff(taus) <= 0))
        tau_duplicates = int(len(taus) - len(np.unique(taus)))
        tau_reach = np.hypot(x[taus] - x[starts], y[taus] - y[starts]) / 1000.0
        if len(tau_reach) != len(frozen) or np.any(tau_reach <= 0) or np.any(taus >= len(g)):
            raise RuntimeError(f"invalid tau confirmation reach: {key}")
        rho_fit = fit_step_distributions(frozen)
        tau_fit = fit_step_distributions(tau_reach)
        rows.append({
            "segment_key": key, "track_id": str(oldrow.track_id), "mean_logchl": float(oldrow.mean_logchl),
            "initial": int(starts[0]), "x": x, "y": y, "candidates": candidates,
            "rho_lengths": frozen, "tau_reach_lengths": tau_reach,
            "rho_support": float(rho_fit["aic_exp"] - rho_fit["aic_lomax"]),
            "tau_reach_support": float(tau_fit["aic_exp"] - tau_fit["aic_lomax"]),
            "rho_median_log_length": float(np.median(np.log(frozen))),
            "tau_reach_median_log_length": float(np.median(np.log(tau_reach))),
        })
        audits.append({
            "segment_key": key, "track_id": str(oldrow.track_id), "events": len(ev),
            "old_n_steps": int(oldrow.n_steps), "event_count_match": len(ev) == int(oldrow.n_steps),
            "rho_coordinate_max_abs_error_km": rho_error,
            "tau_nonincreasing_adjacent_pairs": tau_nonincreasing, "tau_duplicate_fixes": tau_duplicates,
            "tau_confirmation_reach_positive": bool(np.all(tau_reach > 0)),
            "rho_support_old": float(oldrow.delta_aic_lomax), "rho_support_rebuilt": rows[-1]["rho_support"],
            "rho_support_abs_error": abs(float(oldrow.delta_aic_lomax) - rows[-1]["rho_support"]),
            "eligible_records_min": min(map(len, candidates)), "eligible_records_median": float(np.median(list(map(len, candidates)))),
        })
        if done % max(1, len(old) // 20) == 0 or done == len(old):
            progress.emit("build_segments", done, len(old))
    return rows, pd.DataFrame(audits), old


def segment_frame(segments: list[dict[str, Any]], scheme: str) -> pd.DataFrame:
    return pd.DataFrame({
        "segment_key": [s["segment_key"] for s in segments],
        "track_id": [s["track_id"] for s in segments],
        "mean_logchl": [s["mean_logchl"] for s in segments],
        "support": [s[f"{scheme}_support"] for s in segments],
        "median_log_length": [s[f"{scheme}_median_log_length"] for s in segments],
    })


def correlations(frame: pd.DataFrame) -> tuple[float, float]:
    return (
        float(spearmanr(frame.mean_logchl, frame.support).statistic),
        float(spearmanr(frame.mean_logchl, frame.median_log_length).statistic),
    )


def bootstrap(frame: pd.DataFrame, reps: int, seed: int) -> dict[str, Any]:
    tracks = frame.track_id.drop_duplicates().to_numpy(object)
    by = {t: frame[frame.track_id == t] for t in tracks}
    rng = np.random.default_rng(seed)
    observed = correlations(frame)
    values = np.empty((reps, 2), float)
    for i in range(reps):
        sampled = rng.choice(tracks, size=len(tracks), replace=True)
        d = pd.concat([by[t] for t in sampled], ignore_index=True)
        values[i] = correlations(d)
    lo, hi = np.quantile(values, [0.025, 0.975], axis=0)
    return {
        "rho_chl_support": observed[0], "rho_chl_support_ci_low": lo[0], "rho_chl_support_ci_high": hi[0],
        "rho_chl_median_log_length": observed[1], "rho_chl_median_log_length_ci_low": lo[1],
        "rho_chl_median_log_length_ci_high": hi[1], "bootstrap_reps": reps,
    }


def record_replicate(rep: int) -> dict[str, Any]:
    rng = np.random.default_rng(SEED + 10000 + rep)
    rows = []
    for s in _SEGMENTS:
        chosen = np.asarray([c[int(rng.integers(0, len(c)))] for c in s["candidates"]], int)
        lengths = lengths_from_boundaries(s["x"], s["y"], s["initial"], chosen)
        if len(lengths) != len(s["candidates"]) or np.any(lengths <= 0):
            raise RuntimeError(f"invalid record partition at replicate {rep}: {s['segment_key']}")
        fit = fit_step_distributions(lengths)
        rows.append((s["mean_logchl"], fit["aic_exp"] - fit["aic_lomax"], float(np.median(np.log(lengths)))))
    a = np.asarray(rows, float)
    return {
        "replicate": rep,
        "rho_chl_support": float(spearmanr(a[:, 0], a[:, 1]).statistic),
        "rho_chl_median_log_length": float(spearmanr(a[:, 0], a[:, 2]).statistic),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--record-reps", type=int, default=999)
    ap.add_argument("--bootstrap-reps", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=16)
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
        args.record_reps = min(args.record_reps, 19)
        args.bootstrap_reps = min(args.bootstrap_reps, 200)
        args.workers = min(args.workers, 2)
    self_tests()
    progress.emit("start", 0, 1, affinity=affinity, workers=args.workers, smoke=args.smoke)
    global _SEGMENTS
    _SEGMENTS, audit, old = build_segments(args.smoke, progress)
    max_coord = float(audit.rho_coordinate_max_abs_error_km.max())
    max_support = float(audit.rho_support_abs_error.max())
    if not args.smoke and (max_coord >= .01 or max_support >= 1e-6 or not audit.event_count_match.all()):
        raise RuntimeError(f"rho anchor failed coord={max_coord} support={max_support} count={audit.event_count_match.all()}")
    rho_frame = segment_frame(_SEGMENTS, "rho")
    tau_frame = segment_frame(_SEGMENTS, "tau_reach")
    rho_boot = bootstrap(rho_frame, args.bootstrap_reps, SEED + 1)
    tau_boot = bootstrap(tau_frame, args.bootstrap_reps, SEED + 2)
    progress.emit("rho_tau_bootstrap", 2, 2)

    results: list[dict[str, Any]] = []
    if args.workers <= 1:
        for rep in range(args.record_reps):
            results.append(record_replicate(rep))
            if (rep + 1) % max(1, args.record_reps // 20) == 0 or rep + 1 == args.record_reps:
                progress.emit("record_counterfactual", rep + 1, args.record_reps)
    else:
        with futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            jobs = {pool.submit(record_replicate, rep): rep for rep in range(args.record_reps)}
            for done, job in enumerate(futures.as_completed(jobs), 1):
                results.append(job.result())
                if done % max(1, args.record_reps // 20) == 0 or done == args.record_reps:
                    progress.emit("record_counterfactual", done, args.record_reps)
    record = pd.DataFrame(results).sort_values("replicate")
    summary_record = {}
    for metric in ("rho_chl_support", "rho_chl_median_log_length"):
        q = np.quantile(record[metric], [0.025, 0.5, 0.975])
        observed = rho_boot[metric]
        summary_record[metric] = {
            "q025": q[0], "median": q[1], "q975": q[2],
            "rho_observed": observed,
            "one_sided_p_rho_more_negative": float((1 + np.sum(record[metric] <= observed)) / (len(record) + 1)),
        }

    rho_primary = rho_boot["rho_chl_support_ci_high"] < 0 and rho_boot["rho_chl_median_log_length_ci_high"] < 0
    record_both = all(rho_boot[m] < summary_record[m]["median"] for m in summary_record)
    record_tail = any(summary_record[m]["one_sided_p_rho_more_negative"] <= .025 for m in summary_record)
    anchors = max_coord < .01 and max_support < 1e-6 and bool(audit.event_count_match.all()) and bool(audit.tau_confirmation_reach_positive.all())
    passed = bool(rho_primary and record_both and record_tail and anchors)

    audit.to_csv(args.output / "segment_boundary_audit.csv", index=False)
    rho_frame.assign(scheme="rho").to_csv(args.output / "rho_segment_metrics.csv", index=False)
    tau_frame.assign(scheme="tau_overlapping_confirmation_reach").to_csv(args.output / "tau_confirmation_reach_diagnostic.csv", index=False)
    record.to_csv(args.output / "eligible_record_counterfactual_distribution.csv.gz", index=False, compression="gzip")
    payload = {
        "status": "SMOKE_COMPLETE" if args.smoke else "FORMAL_COMPLETE",
        "verdict": "MACRO_TAIL_LENGTH_PARTLY_BOUNDARY_PLACEMENT_SUPPORTED" if passed else "MACRO_BOUNDARY_PLACEMENT_UPGRADE_GATE_FAILED",
        "gate_pass": passed, "rho": rho_boot,
        "tau_overlapping_confirmation_reach_diagnostic": tau_boot,
        "eligible_record_counterfactual": summary_record,
        "gates": {"rho_primary": rho_primary, "record_both_more_negative": record_both,
                  "record_at_least_one_p_le_025": record_tail, "anchors": anchors},
        "invariant": "Event count is identical under rho and eligible-record boundaries; renewal density is therefore unchanged by construction and was not tested as a boundary-placement consequence. Tau marks overlap/cross and are not a renewal partition.",
        "anchor_max_coordinate_error_km": max_coord, "anchor_max_support_error": max_support,
        "segments": len(_SEGMENTS), "tracks": int(old.track_id.nunique()), "record_reps": args.record_reps,
        "workers": args.workers, "affinity": affinity, "cpu1_excluded": 1 not in affinity,
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (CANON, DENSE, OLD, PREREG)},
        "script_sha256": sha256(Path(__file__)),
        "interpretation_limit": "Statistical boundary-placement counterfactual; not an online animal rule, causal intervention, or test of renewal density.",
    }
    (args.output / "final_summary.json").write_text(json.dumps(finite(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("complete", 1, 1, verdict=payload["verdict"])


if __name__ == "__main__":
    main()
