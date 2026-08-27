#!/usr/bin/env python3
"""D104: decompose retrospective CHL endpoint excess under last-passage theory.

This script never reruns RD. It reconstructs strict radial-record fixes inside
frozen RD events, decomposes endpoint excess exactly, and contrasts the frozen
last maximum rho with the frozen drawdown trigger tau.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()

from paper2_core import xy_from_latlon  # noqa: E402


PREREG = ROOT / "metadata/last_record_decomposition_specification_cn.md"
INPUTS = ROOT / "external_inputs/last_record_decomposition"
GOTO_EVENTS = INPUTS / "goto_multiscale_primary_events.csv.gz"
GOTO_CANON = INPUTS / "goto_canonical_events.csv"
GOTO_DENSE = INPUTS / "goto_continuous_paths_with_chl.csv.gz"
RFBO_GPS = INPUTS / "booby_gps_with_daily_chl.csv.gz"
RFBO_RD = INPUTS / "booby_movement_endpoints_all_scales.csv.gz"
RFBO_OLD = ROOT / "data/audit_inputs/rfbo_event_two_tail_metrics.csv.gz"
DEPLOYMENTS = INPUTS / "booby_deployments.csv"
DEFAULT_OUT = ROOT / "results/last_record_decomposition"

RFBO_SCALES = (250, 500, 1000, 2000)
TAILS = ("high", "low", "union")
A_NAMES = tuple(f"{part}_{tail}" for part in ("E", "R", "L") for tail in TAILS)
B_NAMES = tuple(f"{part}_{tail}" for part in ("rho", "tau", "C") for tail in TAILS)
SEED = 1040825
CELL_RES = 1.0 / 24.0


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def holm(values: list[float]) -> list[float]:
    raw = np.asarray(values, float)
    order = np.argsort(np.where(np.isfinite(raw), raw, 1.0))
    out = np.ones(len(raw), float)
    running = 0.0
    for rank, idx in enumerate(order):
        p = raw[idx] if np.isfinite(raw[idx]) else 1.0
        running = max(running, min(1.0, (len(raw) - rank) * p))
        out[idx] = running
    return out.tolist()


class Progress:
    def __init__(self, output: Path):
        self.path = output / "progress.jsonl"
        self.started = time.monotonic()

    def emit(self, stage: str, done: int, total: int, **extra: Any) -> None:
        elapsed = max(time.monotonic() - self.started, 1e-9)
        rate = done / elapsed if done else 0.0
        row = {
            "time": pd.Timestamp.now(tz="America/Los_Angeles").isoformat(),
            "stage": stage, "completed": int(done), "total": int(total),
            "percent": 100.0 * done / max(total, 1), "elapsed_s": elapsed,
            "throughput_per_s": rate,
            "eta_s": (total - done) / rate if rate > 0 and done < total else 0.0,
            **extra,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(finite(row), ensure_ascii=False) + "\n")
        print(
            f"PROGRESS {stage} {done}/{total} ({row['percent']:.1f}%) "
            f"elapsed={elapsed:.1f}s eta={row['eta_s']:.1f}s",
            flush=True,
        )


def add_cell_id(frame: pd.DataFrame, time_col: str, lat_col: str, lon_col: str) -> pd.Series:
    timestamp = pd.to_datetime(frame[time_col], errors="coerce", utc=True, format="mixed")
    lat_bin = np.floor((pd.to_numeric(frame[lat_col], errors="coerce") + 90.0) / CELL_RES).astype("Int64")
    lon_bin = np.floor(((pd.to_numeric(frame[lon_col], errors="coerce") + 180.0) % 360.0) / CELL_RES).astype("Int64")
    return (
        timestamp.dt.strftime("%Y-%m-%d").fillna("NA")
        + "|lat=" + lat_bin.astype(str) + "|lon=" + lon_bin.astype(str)
    )


def tail_matrix(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    ranks = rankdata(values, method="average") / len(values)
    high = ranks > 0.8
    low = ranks <= 0.2
    return np.column_stack([high, low, high | low]).astype(float)


def metrics_a(values: np.ndarray, record_runs: np.ndarray) -> np.ndarray:
    flags = tail_matrix(values)
    overall = flags.mean(axis=0)
    record = flags[np.asarray(record_runs, int)].mean(axis=0)
    endpoint = flags[-1]
    output = np.r_[endpoint - overall, record - overall, endpoint - record]
    if not np.allclose(output[:3], output[3:6] + output[6:9], rtol=0, atol=1e-12):
        raise AssertionError("Family A identity E=R+L failed")
    return output


def metrics_b(values: np.ndarray, rho_run: int) -> np.ndarray:
    flags = tail_matrix(values)
    overall = flags.mean(axis=0)
    rho = flags[int(rho_run)] - overall
    tau = flags[-1] - overall
    contrast = flags[int(rho_run)] - flags[-1]
    output = np.r_[rho, tau, contrast]
    if not np.allclose(output[6:9], output[:3] - output[3:6], rtol=0, atol=1e-12):
        raise AssertionError("Family B identity C=rho-tau failed")
    return output


def strict_record_indices(x: np.ndarray, y: np.ndarray, start: int, endpoint: int) -> np.ndarray:
    origin = np.array([x[start], y[start]], float)
    indices: list[int] = []
    maximum = 0.0
    for index in range(start + 1, endpoint + 1):
        radius = float(np.linalg.norm(np.array([x[index], y[index]]) - origin))
        if radius > maximum:
            maximum = radius
            indices.append(index)
    return np.asarray(indices, dtype=int)


def make_sequence(
    positions: np.ndarray,
    values: np.ndarray,
    cells: np.ndarray,
    start: int,
    end: int,
    marks: dict[str, int],
) -> tuple[np.ndarray, dict[str, int], str]:
    mask = (positions >= start) & (positions <= end) & np.isfinite(values)
    p = positions[mask]
    v = values[mask]
    c = cells[mask]
    if len(v) == 0:
        return np.asarray([], float), {}, "empty_or_no_finite_chl"
    run = np.cumsum(np.r_[True, c[1:] != c[:-1]]) - 1
    mapped: dict[str, int] = {}
    for name, position in marks.items():
        hits = np.flatnonzero(np.isclose(p, float(position), rtol=0, atol=1e-9))
        if len(hits) == 0:
            return np.asarray([], float), {}, f"mark_missing_{name}"
        mapped[name] = int(run[hits[-1]])
    output: list[float] = []
    end_run = mapped.get("end")
    for run_id in range(int(run[-1]) + 1):
        use = run == run_id
        if end_run == run_id:
            exact = use & np.isclose(p, float(end), rtol=0, atol=1e-9)
            output.append(float(v[exact][-1]))
        else:
            output.append(float(np.median(v[use])))
    return np.asarray(output, float), mapped, "ok"


def eligible_values(values: np.ndarray) -> bool:
    return bool(len(values) >= 3 and np.isfinite(values).all() and np.max(values) > np.min(values))


def goto_records(progress: Progress, smoke: bool) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    formal = pd.read_csv(GOTO_EVENTS, low_memory=False)
    formal = formal.loc[formal.scale_m.eq(100) & ~formal.is_terminal.astype(bool)].copy()
    canon_cols = ["event_key", "trigger_orig_idx", "pretrigger_resampled_idx"]
    canon = pd.read_csv(GOTO_CANON, usecols=canon_cols, low_memory=False)
    formal = formal.merge(canon, on="event_key", how="left", validate="one_to_one")
    dense = pd.read_csv(GOTO_DENSE, low_memory=False)
    dense["position"] = np.where(
        dense.orig_idx.notna(), dense.orig_idx,
        dense.edge_start_orig_idx + dense.interpolation_fraction_on_edge
        * (dense.edge_end_orig_idx - dense.edge_start_orig_idx),
    )
    dense["cell_id"] = add_cell_id(dense, "Time", "Lat", "Lon")
    keys = formal[["track_id", "segment_id"]].drop_duplicates()
    if smoke:
        keys = keys.head(2)
        keep = set(map(tuple, keys.itertuples(index=False, name=None)))
        formal = formal[[tuple(x) in keep for x in formal[["track_id", "segment_id"]].itertuples(index=False, name=None)]]
    groups = {(str(a), int(b)): g.sort_values("position", kind="stable") for (a, b), g in dense.groupby(["track_id", "segment_id"], sort=False)}
    a_rows: list[dict[str, Any]] = []
    b_rows: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    grouped_events = list(formal.groupby(["track_id", "segment_id"], sort=False))
    for number, ((track, segment), events) in enumerate(grouped_events, 1):
        g = groups[(str(track), int(segment))]
        originals = g.loc[g.orig_idx.notna()].sort_values("orig_idx", kind="stable")
        original_index = originals.orig_idx.to_numpy(int)
        if not np.array_equal(original_index, np.arange(len(originals))):
            raise RuntimeError(f"Goto original index is not contiguous: {track}/{segment}")
        x, y = xy_from_latlon(originals.Lat.to_numpy(float), originals.Lon.to_numpy(float))
        pos = g.position.to_numpy(float)
        val = g.logCHL.to_numpy(float)
        cell = g.cell_id.astype(str).to_numpy(object)
        for event in events.sort_values("step_id").itertuples(index=False):
            start = int(event.start_index); rho = int(event.endpoint_index); tau = int(event.trigger_orig_idx)
            rec = strict_record_indices(x, y, start, rho)
            geometry_ok = bool(len(rec) and rec[-1] == rho)
            values_a, marks_a, status_a = make_sequence(
                pos, val, cell, start, rho,
                {"end": rho, **{f"record_{i}": int(index) for i, index in enumerate(rec)}},
            )
            record_runs = np.unique([marks_a[key] for key in marks_a if key.startswith("record_")]) if status_a == "ok" else np.asarray([], int)
            a_eligible = geometry_ok and status_a == "ok" and eligible_values(values_a)
            a_decomp = a_eligible and len(record_runs) >= 2
            if a_decomp:
                a_rows.append({
                    "dataset": "goto", "scale_m": 100, "family": "A", "event_id": str(event.event_key),
                    "segment_key": f"{track}|{int(segment)}", "cluster": str(event.cluster_id),
                    "values": values_a, "position": record_runs.astype(int),
                    "observed": metrics_a(values_a, record_runs), "rho_tau_distinct_cell": None,
                })
            values_b, marks_b, status_b = make_sequence(
                pos, val, cell, start, tau, {"rho": rho, "end": tau},
            )
            b_eligible = status_b == "ok" and eligible_values(values_b)
            if b_eligible:
                rho_run = int(marks_b["rho"])
                b_rows.append({
                    "dataset": "goto", "scale_m": 100, "family": "B", "event_id": str(event.event_key),
                    "segment_key": f"{track}|{int(segment)}", "cluster": str(event.cluster_id),
                    "values": values_b, "position": rho_run,
                    "observed": metrics_b(values_b, rho_run),
                    "rho_tau_distinct_cell": bool(rho_run != len(values_b) - 1),
                })
            audit.append({
                "dataset": "goto", "scale_m": 100, "event_id": str(event.event_key),
                "segment_key": f"{track}|{int(segment)}", "cluster": str(event.cluster_id),
                "geometry_endpoint_is_last_record": geometry_ok, "a_status": status_a,
                "a_runs": len(values_a), "a_record_runs": len(record_runs), "a_eligible": a_eligible,
                "a_decomposition_eligible": a_decomp, "b_status": status_b, "b_runs": len(values_b),
                "b_eligible": b_eligible,
                "rho_tau_distinct_cell": bool(status_b == "ok" and marks_b.get("rho") != marks_b.get("end")),
            })
        progress.emit("goto_reconstruct", number, len(grouped_events), a_events=len(a_rows), b_events=len(b_rows))
    return a_rows, b_rows, audit


def rfbo_records(progress: Progress, smoke: bool) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    gps = pd.read_csv(RFBO_GPS, low_memory=False)
    gps["cell_id"] = add_cell_id(gps, "Timestamp", "lat", "lon")
    rd = pd.read_csv(RFBO_RD, low_memory=False)
    rd = rd.loc[rd.delta_m.isin(RFBO_SCALES)].copy()
    deployments = pd.read_csv(DEPLOYMENTS, usecols=["DeployID", "BandNum"])
    band_map = deployments.set_index("DeployID").BandNum.astype(str).to_dict()
    frozen_d28 = pd.read_csv(RFBO_OLD, usecols=["event_id", "delta_m"])
    frozen_d28_keys = set(
        (int(scale), str(event_id))
        for scale, event_id in frozen_d28[["delta_m", "event_id"]].itertuples(index=False, name=None)
    )
    if smoke:
        keep_segments = rd.segment_id.drop_duplicates().head(2)
        rd = rd.loc[rd.segment_id.isin(keep_segments)]
    groups = {str(key): g.sort_values("segment_pos", kind="stable").reset_index(drop=True) for key, g in gps.groupby("segment_id", sort=False)}
    a_rows: list[dict[str, Any]] = []
    b_rows: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    grouped_events = list(rd.groupby(["delta_m", "segment_id"], sort=False))
    for number, ((scale, segment), events) in enumerate(grouped_events, 1):
        g = groups[str(segment)]
        positions = g.segment_pos.to_numpy(int)
        if not np.array_equal(positions, np.arange(len(g))):
            raise RuntimeError(f"RFBO segment_pos is not contiguous: {segment}")
        lat0 = float(g.lat.mean())
        x = (g.lon.to_numpy(float) - float(g.lon.iloc[0])) * 111.32 * math.cos(math.radians(lat0))
        y = (g.lat.to_numpy(float) - float(g.lat.iloc[0])) * 110.574
        values = g.log_chl.to_numpy(float)
        cells = g.cell_id.astype(str).to_numpy(object)
        cluster = str(band_map.get(int(g.DeployID.iloc[0]), f"DEPLOY_{int(g.DeployID.iloc[0])}"))
        for event in events.sort_values("start_idx").itertuples(index=False):
            start = int(event.start_idx); rho = int(event.endpoint_idx); tau = int(event.trigger_idx)
            rec = strict_record_indices(x, y, start, rho)
            geometry_ok = bool(len(rec) and rec[-1] == rho)
            values_a, marks_a, status_a = make_sequence(
                positions, values, cells, start, rho,
                {"end": rho, **{f"record_{i}": int(index) for i, index in enumerate(rec)}},
            )
            record_runs = np.unique([marks_a[key] for key in marks_a if key.startswith("record_")]) if status_a == "ok" else np.asarray([], int)
            a_eligible = geometry_ok and status_a == "ok" and eligible_values(values_a)
            in_frozen_d28 = (int(scale), str(event.event_id)) in frozen_d28_keys
            a_decomp = a_eligible and len(record_runs) >= 2 and in_frozen_d28
            common = {
                "dataset": "rfbo", "scale_m": int(scale), "event_id": str(event.event_id),
                "segment_key": str(segment), "cluster": cluster,
            }
            if a_decomp:
                a_rows.append({**common, "family": "A", "values": values_a,
                               "position": record_runs.astype(int), "observed": metrics_a(values_a, record_runs),
                               "rho_tau_distinct_cell": None})
            values_b, marks_b, status_b = make_sequence(
                positions, values, cells, start, tau, {"rho": rho, "end": tau},
            )
            b_eligible = status_b == "ok" and eligible_values(values_b)
            if b_eligible:
                rho_run = int(marks_b["rho"])
                b_rows.append({**common, "family": "B", "values": values_b, "position": rho_run,
                               "observed": metrics_b(values_b, rho_run),
                               "rho_tau_distinct_cell": bool(rho_run != len(values_b) - 1)})
            audit.append({
                **common, "geometry_endpoint_is_last_record": geometry_ok, "a_status": status_a,
                "a_runs": len(values_a), "a_record_runs": len(record_runs), "a_eligible": a_eligible,
                "a_in_frozen_d28_universe": in_frozen_d28,
                "a_decomposition_eligible": a_decomp, "b_status": status_b, "b_runs": len(values_b),
                "b_eligible": b_eligible,
                "rho_tau_distinct_cell": bool(status_b == "ok" and marks_b.get("rho") != marks_b.get("end")),
            })
        progress.emit("rfbo_reconstruct", number, len(grouped_events), a_events=len(a_rows), b_events=len(b_rows))
    return a_rows, b_rows, audit


def shifted_segment_sums(
    records: list[dict[str, Any]], family: str, reps: int, rng: np.random.Generator,
) -> tuple[np.ndarray, int, int]:
    master: list[float] = []
    descriptors: list[tuple[int, int, Any]] = []
    for record in records:
        values = np.asarray(record["values"], float)
        master.extend(values.tolist())
        descriptors.append((len(master) - 1, len(values), record["position"]))
    array = np.asarray(master, float)
    n_metrics = len(A_NAMES if family == "A" else B_NAMES)
    if len(array) < 2:
        raise RuntimeError("segment master has fewer than two CHL tokens")
    metric: Callable[[np.ndarray, Any], np.ndarray] = metrics_a if family == "A" else metrics_b
    cache: dict[int, np.ndarray | None] = {}
    outputs = np.zeros((reps, n_metrics), float)
    attempts = 0
    for replicate in range(reps):
        while True:
            attempts += 1
            if attempts > reps * max(1000, 20 * len(array)):
                raise RuntimeError("could not sample a valid nonzero common phase")
            offset = int(rng.integers(1, len(array)))
            if offset not in cache:
                total = np.zeros(n_metrics, float)
                valid = True
                for end, width, position in descriptors:
                    pseudo_end = (end - offset) % len(array)
                    index = (pseudo_end - np.arange(width - 1, -1, -1, dtype=int)) % len(array)
                    window = array[index]
                    if np.max(window) <= np.min(window):
                        valid = False
                        break
                    total += metric(window, position)
                cache[offset] = total if valid else None
            value = cache[offset]
            if value is not None:
                outputs[replicate] = value
                break
    return outputs, len(cache), int(sum(value is not None for value in cache.values()))


def summarize_family(
    records: list[dict[str, Any]], family: str, reps: int, bootstrap_reps: int,
    rng: np.random.Generator, progress: Progress,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    names = A_NAMES if family == "A" else B_NAMES
    output_rows: list[dict[str, Any]] = []
    null_audit: list[dict[str, Any]] = []
    groups = list(pd.DataFrame([
        {"dataset": r["dataset"], "scale_m": r["scale_m"]} for r in records
    ]).drop_duplicates().itertuples(index=False))
    for group_number, group in enumerate(groups, 1):
        use = [r for r in records if r["dataset"] == group.dataset and r["scale_m"] == group.scale_m]
        segment_map: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for record in use:
            segment_map[str(record["segment_key"])].append(record)
        cluster_observed_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(len(names), float))
        cluster_count: dict[str, int] = defaultdict(int)
        cluster_null_sum: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((reps, len(names)), float))
        for segment_number, (segment, segment_records) in enumerate(segment_map.items(), 1):
            clusters = {str(record["cluster"]) for record in segment_records}
            if len(clusters) != 1:
                raise RuntimeError(f"segment maps to multiple clusters: {segment}")
            cluster = next(iter(clusters))
            observed_sum = np.sum([np.asarray(record["observed"], float) for record in segment_records], axis=0)
            shifted, attempted, valid = shifted_segment_sums(segment_records, family, reps, rng)
            cluster_observed_sum[cluster] += observed_sum
            cluster_null_sum[cluster] += shifted
            cluster_count[cluster] += len(segment_records)
            null_audit.append({
                "dataset": group.dataset, "scale_m": int(group.scale_m), "family": family,
                "segment_key": segment, "cluster": cluster, "events": len(segment_records),
                "master_tokens": int(sum(len(record["values"]) for record in segment_records)),
                "attempted_unique_offsets": attempted, "valid_unique_offsets": valid,
            })
            if segment_number % max(1, len(segment_map) // 10) == 0 or segment_number == len(segment_map):
                progress.emit(
                    f"null_{family}_{group.dataset}_{int(group.scale_m)}",
                    segment_number, len(segment_map), events=len(use),
                )
        clusters = sorted(cluster_count)
        unit = np.row_stack([cluster_observed_sum[c] / cluster_count[c] for c in clusters])
        observed = unit.mean(axis=0)
        null = np.mean(np.stack([cluster_null_sum[c] / cluster_count[c] for c in clusters], axis=1), axis=1)
        chosen = rng.integers(0, len(unit), size=(bootstrap_reps, len(unit)))
        bootstrap = unit[chosen].mean(axis=1)
        lo, hi = np.quantile(bootstrap, [0.025, 0.975], axis=0)
        distinct = [bool(r["rho_tau_distinct_cell"]) for r in use if r["rho_tau_distinct_cell"] is not None]
        for index, name in enumerate(names):
            output_rows.append({
                "dataset": group.dataset, "scale_m": int(group.scale_m), "family": family,
                "estimand": name, "events": len(use), "units": len(clusters), "segments": len(segment_map),
                "observed_unit_equal": float(observed[index]),
                "bootstrap_ci_low": float(lo[index]), "bootstrap_ci_high": float(hi[index]),
                "phase_null_mean": float(null[:, index].mean()),
                "phase_p_one_sided": float((1 + np.sum(null[:, index] >= observed[index])) / (reps + 1)),
                "phase_reps": reps, "bootstrap_reps": bootstrap_reps,
                "rho_tau_distinct_cell_fraction": float(np.mean(distinct)) if distinct else np.nan,
            })
        progress.emit("summarize_family", group_number, len(groups), family=family, dataset=group.dataset, scale_m=int(group.scale_m))
    summary = pd.DataFrame(output_rows)
    summary["holm_phase_p"] = np.nan
    for (dataset, fam), idx in summary.groupby(["dataset", "family"]).groups.items():
        if fam == "A":
            tested = [i for i in idx if str(summary.loc[i, "estimand"]).startswith(("R_", "L_"))]
        else:
            tested = list(idx)
        summary.loc[tested, "holm_phase_p"] = holm(summary.loc[tested, "phase_p_one_sided"].tolist())
    summary["positive_support"] = (
        summary.bootstrap_ci_low.gt(0) & summary.holm_phase_p.le(0.05)
    )
    return summary, pd.DataFrame(null_audit)


def event_frame(records: list[dict[str, Any]], names: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for record in records:
        row = {key: value for key, value in record.items() if key not in {"values", "position", "observed"}}
        row["n_runs"] = len(record["values"])
        if record["family"] == "A":
            row["n_record_runs"] = len(record["position"])
        else:
            row["rho_run"] = int(record["position"])
        for name, value in zip(names, record["observed"]):
            row[name] = float(value)
        rows.append(row)
    return pd.DataFrame(rows)


def anchor_audit(a_events: pd.DataFrame) -> dict[str, Any]:
    goto_old = pd.read_csv(GOTO_EVENTS)
    goto_old = goto_old.loc[goto_old.scale_m.eq(100), ["event_id", "high20_excess", "low20_excess", "extreme20_excess"]]
    goto = a_events.loc[a_events.dataset.eq("goto")].merge(goto_old, on="event_id", how="left", validate="one_to_one")
    rfbo_old = pd.read_csv(RFBO_OLD, usecols=["event_id", "delta_m", "high20_excess", "low20_excess", "extreme20_excess"])
    rfbo = a_events.loc[a_events.dataset.eq("rfbo")].merge(
        rfbo_old, left_on=["event_id", "scale_m"], right_on=["event_id", "delta_m"],
        how="left", validate="one_to_one",
    )
    result: dict[str, Any] = {}
    for dataset, frame in (("goto", goto), ("rfbo", rfbo)):
        diffs = []
        for new, old in (("E_high", "high20_excess"), ("E_low", "low20_excess"), ("E_union", "extreme20_excess")):
            diff = np.abs(frame[new].to_numpy(float) - frame[old].to_numpy(float))
            diffs.extend(diff[np.isfinite(diff)].tolist())
        result[dataset] = {
            "rows": len(frame), "old_rows_matched": int(frame.high20_excess.notna().sum()),
            "maximum_absolute_E_difference": float(max(diffs)) if diffs else None,
            "exact_within_1e_12_fraction": float(np.mean(np.asarray(diffs) <= 1e-12)) if diffs else None,
        }
    return result


def verdict(summary: pd.DataFrame) -> dict[str, Any]:
    answer: dict[str, Any] = {"family_a": {}, "family_b": {}}
    for tail in TAILS:
        component_pass: dict[str, bool] = {}
        for part in ("R", "L"):
            name = f"{part}_{tail}"
            goto = summary[(summary.dataset == "goto") & (summary.family == "A") & (summary.estimand == name)]
            rfbo = summary[(summary.dataset == "rfbo") & (summary.family == "A") & (summary.estimand == name)]
            scales = sorted(rfbo.loc[rfbo.positive_support, "scale_m"].astype(int).tolist())
            adjacent = any(a in scales and b in scales for a, b in zip(RFBO_SCALES[:-1], RFBO_SCALES[1:]))
            component_pass[part] = bool(len(goto) == 1 and bool(goto.positive_support.iloc[0]) and adjacent)
        if component_pass["R"] and component_pass["L"]:
            value = "BOTH_COMPONENTS_REPLICATED"
        elif component_pass["R"]:
            value = "RECORD_SAMPLING_DOMINANT"
        elif component_pass["L"]:
            value = "LAST_RECORD_SELECTION_DOMINANT"
        else:
            value = "NO_REPLICATED_COMPONENT_LOCALIZATION"
        answer["family_a"][tail] = {"verdict": value, "component_cross_system_pass": component_pass}
        name = f"C_{tail}"
        goto = summary[(summary.dataset == "goto") & (summary.family == "B") & (summary.estimand == name)]
        rfbo = summary[(summary.dataset == "rfbo") & (summary.family == "B") & (summary.estimand == name)]
        scales = sorted(rfbo.loc[rfbo.positive_support, "scale_m"].astype(int).tolist())
        adjacent = any(a in scales and b in scales for a, b in zip(RFBO_SCALES[:-1], RFBO_SCALES[1:]))
        replicated = bool(len(goto) == 1 and bool(goto.positive_support.iloc[0]) and adjacent)
        answer["family_b"][tail] = {
            "verdict": "LAST_MAXIMUM_LOCALIZATION_REPLICATED" if replicated else "NO_REPLICATED_LAST_MAXIMUM_LOCALIZATION",
            "goto_pass": bool(len(goto) == 1 and bool(goto.positive_support.iloc[0])),
            "rfbo_pass_scales": scales,
        }
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--prereg", type=Path, default=PREREG)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--phase-reps", type=int, default=199)
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.output)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 is forbidden; affinity={affinity}")
    if args.smoke:
        args.phase_reps = min(args.phase_reps, 19)
        args.bootstrap_reps = min(args.bootstrap_reps, 200)
    progress.emit("start", 0, 1, affinity=affinity, smoke=args.smoke)
    rng = np.random.default_rng(SEED)
    goto_a, goto_b, goto_audit = goto_records(progress, args.smoke)
    rfbo_a, rfbo_b, rfbo_audit = rfbo_records(progress, args.smoke)
    a_records = goto_a + rfbo_a
    b_records = goto_b + rfbo_b
    if not a_records or not b_records:
        raise RuntimeError("no eligible D104 records")
    a_frame = event_frame(a_records, A_NAMES)
    b_frame = event_frame(b_records, B_NAMES)
    max_a_identity = float(np.max(np.abs(a_frame[[f"E_{t}" for t in TAILS]].to_numpy() - a_frame[[f"R_{t}" for t in TAILS]].to_numpy() - a_frame[[f"L_{t}" for t in TAILS]].to_numpy())))
    max_b_identity = float(np.max(np.abs(b_frame[[f"C_{t}" for t in TAILS]].to_numpy() - b_frame[[f"rho_{t}" for t in TAILS]].to_numpy() + b_frame[[f"tau_{t}" for t in TAILS]].to_numpy())))
    if max_a_identity > 1e-12 or max_b_identity > 1e-12:
        raise RuntimeError(f"identity audit failed A={max_a_identity} B={max_b_identity}")
    summary_a, null_a = summarize_family(a_records, "A", args.phase_reps, args.bootstrap_reps, rng, progress)
    summary_b, null_b = summarize_family(b_records, "B", args.phase_reps, args.bootstrap_reps, rng, progress)
    summary = pd.concat([summary_a, summary_b], ignore_index=True)
    audit = pd.DataFrame(goto_audit + rfbo_audit)
    coverage = (
        audit.groupby(["dataset", "scale_m"], as_index=False)
        .agg(
            frozen_events=("event_id", "size"), geometry_pass=("geometry_endpoint_is_last_record", "sum"),
            a_chl_eligible=("a_eligible", "sum"), a_decomposition_eligible=("a_decomposition_eligible", "sum"),
            b_chl_eligible=("b_eligible", "sum"), rho_tau_distinct_cell=("rho_tau_distinct_cell", "sum"),
        )
    )
    coverage["rho_tau_distinct_cell_fraction"] = coverage.rho_tau_distinct_cell / coverage.b_chl_eligible
    anchor = anchor_audit(a_frame) if not args.smoke else {"status": "not_run_in_smoke"}
    final_verdict = verdict(summary)
    a_frame.to_csv(args.output / "family_a_event_metrics.csv.gz", index=False, compression="gzip")
    b_frame.to_csv(args.output / "family_b_event_metrics.csv.gz", index=False, compression="gzip")
    audit.to_csv(args.output / "event_structure_audit.csv.gz", index=False, compression="gzip")
    coverage.to_csv(args.output / "coverage.csv", index=False)
    summary.to_csv(args.output / "formal_summary.csv", index=False)
    pd.concat([null_a, null_b], ignore_index=True).to_csv(args.output / "phase_null_audit.csv.gz", index=False, compression="gzip")
    payload = {
        "status": "SMOKE_COMPLETE" if args.smoke else "FORMAL_COMPLETE",
        "affinity": affinity, "cpu1_excluded": 1 not in affinity,
        "phase_reps": args.phase_reps, "bootstrap_reps": args.bootstrap_reps,
        "family_a_identity_max_abs_error": max_a_identity,
        "family_b_identity_max_abs_error": max_b_identity,
        "anchor_reproduction": anchor, "verdict": final_verdict,
        "coverage": coverage.to_dict("records"), "summary": summary.to_dict("records"),
        "input_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in [GOTO_EVENTS, GOTO_CANON, GOTO_DENSE, RFBO_GPS, RFBO_RD, RFBO_OLD, DEPLOYMENTS, args.prereg]},
        "script_sha256": sha256(Path(__file__)),
        "interpretation_limit": "Retrospective last-passage decomposition and rho-vs-tau event localization; not sensory cue, prey truth, online hazard, Levy generation, or causal identification.",
    }
    (args.output / "final_summary.json").write_text(json.dumps(finite(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("complete", 1, 1, status=payload["status"], verdict=final_verdict)


if __name__ == "__main__":
    main()
