#!/usr/bin/env python3
"""D107: short-tailed shearwater last-passage CHL x foraging/rest.

Locked specification:
NODE_STATE/D107_SHORT_TAILED_LAST_PASSAGE_CHL_FORAGING_REST_PREREG_20260825_CN.md

This script uses frozen phase-0 RD-P endpoints. It never runs an RD detector.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
import xarray as xr


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
import shearwater_behavior_bridge as short  # noqa: E402
import last_record_decomposition as last_record  # noqa: E402


PREREG = ROOT / "metadata/shearwater_behavior_specification_cn.md"
PLAN = ROOT / "external_inputs/shearwater_behavior"
CATALOG = PLAN / "geometry_event_catalog.csv.gz"
MANIFEST = PLAN / "copernicus_chl_download_manifest.csv"
ENDPOINTS = ROOT / "results/shearwater_behavior_bridge/behavior_blind_five_minute_rd_endpoints.csv"
DEFAULT_OUT = ROOT / "results/shearwater_behavior"
SCALES = (500, 1000, 2000, 5000)
METRICS = tuple(last_record.A_NAMES)
SEED = 1070825


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
    order = np.argsort(np.where(np.isfinite(raw), raw, 1.0), kind="stable")
    output = np.ones(len(raw), float)
    running = 0.0
    for rank, index in enumerate(order):
        value = raw[index] if np.isfinite(raw[index]) else 1.0
        running = max(running, min(1.0, (len(raw) - rank) * value))
        output[index] = running
    return output.tolist()


class Progress:
    def __init__(self, output: Path) -> None:
        self.output = output
        self.path = output / "progress.jsonl"
        self.started = time.monotonic()

    def emit(self, stage: str, done: int, total: int, **extra: Any) -> None:
        elapsed = max(time.monotonic() - self.started, 1e-9)
        rate = done / elapsed if done else 0.0
        row = {
            "stage": stage, "completed": int(done), "total": int(total),
            "percent": round(100.0 * done / max(total, 1), 3),
            "elapsed_s": round(elapsed, 3), "throughput_per_s": round(rate, 5),
            "eta_s": round((total - done) / rate, 3) if rate > 0 and done < total else 0.0,
            **extra,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(finite(row), ensure_ascii=False) + "\n")
        print(json.dumps(finite(row), ensure_ascii=False), flush=True)


def audit_and_load_chl(progress: Progress) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    manifest = pd.read_csv(MANIFEST)
    expected = [ROOT / str(path) for path in manifest.output_path]
    missing = [str(path) for path in expected if not path.exists()]
    if missing:
        raise RuntimeError(f"missing CHL files ({len(missing)}): {missing[:3]}")

    file_rows: list[dict[str, Any]] = []
    cell_frames: list[pd.DataFrame] = []
    for number, (row, path) in enumerate(zip(manifest.itertuples(index=False), expected), 1):
        with xr.open_dataset(path) as ds:
            if "CHL" not in ds or not {"time", "latitude", "longitude"}.issubset(ds.coords):
                raise RuntimeError(f"invalid CHL schema: {path}")
            dates = pd.to_datetime(ds.time.values).strftime("%Y-%m-%d").tolist()
            if dates != [str(row.start_date)]:
                raise RuntimeError(f"CHL date mismatch: {path}: {dates} != {row.start_date}")
            values = np.asarray(ds.CHL.isel(time=0).values, float)
            lat = np.asarray(ds.latitude.values, float)
            lon = np.asarray(ds.longitude.values, float)
            if values.shape != (len(lat), len(lon)):
                raise RuntimeError(f"CHL shape mismatch: {path}")
            yy, xx = np.meshgrid(lat, lon, indexing="ij")
            cells = pd.DataFrame({
                "date": str(row.start_date), "grid_lat": yy.ravel(), "grid_lon": xx.ravel(),
                "chl": values.ravel(), "source_file": path.name,
            })
            cells = cells.loc[np.isfinite(cells.chl) & (cells.chl > 0)].copy()
            cell_frames.append(cells)
            file_rows.append({
                "file": str(path.relative_to(ROOT)), "sha256": sha256(path), "size_bytes": path.stat().st_size,
                "date": dates[0], "n_lat": len(lat), "n_lon": len(lon),
                "finite_positive_cells": len(cells), "chl_min": float(cells.chl.min()),
                "chl_max": float(cells.chl.max()), "lat_min": float(lat.min()), "lat_max": float(lat.max()),
                "lon_min": float(lon.min()), "lon_max": float(lon.max()),
            })
        progress.emit("chl_file_audit", number, len(expected), file=path.name)

    all_cells = pd.concat(cell_frames, ignore_index=True)
    all_cells["lat_key"] = all_cells.grid_lat.round(6)
    all_cells["lon_key"] = all_cells.grid_lon.round(6)
    overlap = all_cells.groupby(["date", "lat_key", "lon_key"], as_index=False).agg(
        copies=("chl", "size"), value_min=("chl", "min"), value_max=("chl", "max"),
        grid_lat=("grid_lat", "first"), grid_lon=("grid_lon", "first"), chl=("chl", "first"),
    )
    overlap["difference"] = overlap.value_max - overlap.value_min
    max_difference = float(overlap.difference.max()) if len(overlap) else math.nan
    if np.isfinite(max_difference) and max_difference > 1e-6:
        raise RuntimeError(f"overlapping CHL tiles disagree: max difference={max_difference}")
    unique = overlap[["date", "grid_lat", "grid_lon", "chl"]].copy()
    audit = {
        "manifest_rows": len(manifest), "files_found": len(expected),
        "unique_finite_positive_cells": len(unique),
        "overlap_cells": int((overlap.copies > 1).sum()),
        "maximum_overlap_value_difference": max_difference,
    }
    return unique, pd.DataFrame(file_rows), audit


def nearest_chl(points: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    output = points.copy()
    output["date"] = pd.to_datetime(output.timestamp).dt.strftime("%Y-%m-%d")
    output["grid_lat"] = np.nan
    output["grid_lon"] = np.nan
    output["chl"] = np.nan
    output["grid_distance_deg"] = np.nan
    by_date = {date: frame.reset_index(drop=True) for date, frame in cells.groupby("date", sort=False)}
    for date, index in output.groupby("date", sort=False).groups.items():
        if date not in by_date:
            continue
        grid = by_date[date]
        glat = grid.grid_lat.to_numpy(float)
        glon = grid.grid_lon.to_numpy(float)
        for row_index in index:
            lat = float(output.at[row_index, "lat"]); lon = float(output.at[row_index, "lon"])
            scale = math.cos(math.radians(lat))
            dist2 = (glat - lat) ** 2 + ((glon - lon) * scale) ** 2
            nearest = int(np.argmin(dist2))
            output.at[row_index, "grid_lat"] = glat[nearest]
            output.at[row_index, "grid_lon"] = glon[nearest]
            output.at[row_index, "chl"] = float(grid.chl.iloc[nearest])
            output.at[row_index, "grid_distance_deg"] = math.sqrt(float(dist2[nearest]))
    return output


def reconstruct_events(cells: pd.DataFrame, progress: Progress) -> tuple[list[dict[str, Any]], pd.DataFrame, dict[str, Any]]:
    high_runs, _, input_audit = short.load_highres(0)
    supports = short.build_support_runs(high_runs)[0]
    run_map = {str(frame.run_key.iloc[0]): frame.reset_index(drop=True) for frame in supports}
    catalog = pd.read_csv(CATALOG, low_memory=False)
    catalog = catalog.loc[catalog.geometry_eligible.astype(bool) & catalog.scale_m.isin(SCALES)].copy()

    records: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for number, event in enumerate(catalog.itertuples(index=False), 1):
        run = run_map[str(event.run_key)]
        start = int(event.origin_index_support_run); end = int(event.endpoint_index_support_run)
        step = run.iloc[start:end + 1].copy().reset_index(drop=True)
        matched = nearest_chl(step[["timestamp", "lat", "lon"]], cells)
        status = "ok"
        if matched.chl.isna().any():
            status = "missing_chl"
        elif float(matched.grid_distance_deg.max()) > 0.04:
            status = "nearest_grid_too_far"
        if status == "ok":
            cell_id = (
                matched.date.astype(str) + "|" + matched.grid_lat.round(6).astype(str)
                + "|" + matched.grid_lon.round(6).astype(str)
            ).to_numpy(object)
            cell_run = np.cumsum(np.r_[True, cell_id[1:] != cell_id[:-1]]) - 1
            values = []
            for run_id in range(int(cell_run[-1]) + 1):
                values.append(float(np.median(np.log(matched.loc[cell_run == run_id, "chl"].to_numpy(float)))))
            values_array = np.asarray(values, float)
            rec = last_record.strict_record_indices(
                step.x_m.to_numpy(float), step.y_m.to_numpy(float), 0, len(step) - 1
            )
            record_runs = np.unique(cell_run[rec]) if len(rec) else np.asarray([], int)
            if not len(rec) or int(rec[-1]) != len(step) - 1:
                status = "endpoint_not_last_record"
            elif len(values_array) < 3:
                status = "fewer_than_3_native_cell_runs"
            elif not np.isfinite(values_array).all() or np.max(values_array) <= np.min(values_array):
                status = "constant_or_nonfinite_chl"
            elif len(record_runs) < 2:
                status = "fewer_than_2_record_cell_runs"
        if status == "ok":
            observed = last_record.metrics_a(values_array, record_runs)
            row = {
                "dataset": "short_tailed_shearwater", "scale_m": int(event.scale_m),
                "event_id": str(event.event_id), "run_key": str(event.run_key),
                "individual": str(event.individual), "deployment": str(event.deployment),
                "origin_index_support_run": start, "endpoint_index_support_run": end,
                "values": values_array, "position": record_runs.astype(int), "observed": observed,
                "future_specificity": float(run.loc[end, "FUTURE_5MIN__specificity"]),
                "future_forage": float(run.loc[end, "FUTURE_5MIN__forage"]),
                "future_rest": float(run.loc[end, "FUTURE_5MIN__rest"]),
                "prior_specificity": float(run.loc[end, "PRIOR_5MIN__specificity"]),
                "prior_forage": float(run.loc[end, "PRIOR_5MIN__forage"]),
                "prior_rest": float(run.loc[end, "PRIOR_5MIN__rest"]),
                "future_minus_prior_specificity": float(
                    run.loc[end, "FUTURE_5MIN__specificity"] - run.loc[end, "PRIOR_5MIN__specificity"]
                ),
                "behavior_class": "forage_dominant" if event.future_specificity > 0 else (
                    "rest_dominant" if event.future_specificity < 0 else "tie"
                ),
            }
            records.append(row)
            n_cell_runs = len(values_array); n_record_runs = len(record_runs)
        else:
            n_cell_runs = 0; n_record_runs = 0
        audit_rows.append({
            "event_id": str(event.event_id), "scale_m": int(event.scale_m), "run_key": str(event.run_key),
            "individual": str(event.individual), "status": status,
            "support_points": len(step), "native_cell_runs": n_cell_runs,
            "record_cell_runs": n_record_runs,
            "maximum_nearest_grid_distance_deg": float(matched.grid_distance_deg.max()) if matched.grid_distance_deg.notna().any() else math.nan,
        })
        if number % 25 == 0 or number == len(catalog):
            progress.emit("event_reconstruction", number, len(catalog), eligible=len(records))

    audit = pd.DataFrame(audit_rows)
    if not records:
        raise RuntimeError("no eligible D107 events")
    return records, audit, input_audit


def phase_event_metrics(
    records: list[dict[str, Any]], reps: int, rng: np.random.Generator, progress: Progress,
) -> tuple[np.ndarray, pd.DataFrame]:
    output = np.zeros((reps, len(records), len(METRICS)), float)
    groups: dict[tuple[int, str], list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        groups[(int(record["scale_m"]), str(record["run_key"]))].append(index)
    null_rows: list[dict[str, Any]] = []
    for number, ((scale, run_key), event_indices) in enumerate(groups.items(), 1):
        master: list[float] = []
        descriptors: list[tuple[int, int, np.ndarray, int]] = []
        for event_index in event_indices:
            record = records[event_index]
            values = np.asarray(record["values"], float)
            master.extend(values.tolist())
            descriptors.append((len(master) - 1, len(values), np.asarray(record["position"], int), event_index))
        array = np.asarray(master, float)
        if len(array) < 2:
            raise RuntimeError("phase master has fewer than two tokens")
        cache: dict[int, list[np.ndarray] | None] = {}
        attempts = 0
        for replicate in range(reps):
            while True:
                attempts += 1
                if attempts > reps * max(1000, 20 * len(array)):
                    raise RuntimeError(f"could not sample valid phase: {scale}, {run_key}")
                offset = int(rng.integers(1, len(array)))
                if offset not in cache:
                    shifted: list[np.ndarray] = []
                    valid = True
                    for end, width, position, _ in descriptors:
                        pseudo_end = (end - offset) % len(array)
                        index = (pseudo_end - np.arange(width - 1, -1, -1, dtype=int)) % len(array)
                        window = array[index]
                        if np.max(window) <= np.min(window):
                            valid = False
                            break
                        shifted.append(last_record.metrics_a(window, position))
                    cache[offset] = shifted if valid else None
                shifted = cache[offset]
                if shifted is not None:
                    for metric_values, descriptor in zip(shifted, descriptors):
                        output[replicate, descriptor[3], :] = metric_values
                    break
        null_rows.append({
            "scale_m": scale, "run_key": run_key, "events": len(event_indices),
            "master_tokens": len(array), "attempted_draws": attempts,
            "unique_offsets_checked": len(cache),
            "valid_unique_offsets": int(sum(value is not None for value in cache.values())),
        })
        if number % 10 == 0 or number == len(groups):
            progress.emit("common_phase_null", number, len(groups), reps=reps)
    return output, pd.DataFrame(null_rows)


def bootstrap_units(unit: np.ndarray, reps: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if unit.ndim != 2 or len(unit) < 2:
        return np.full(unit.shape[1], np.nan), np.full(unit.shape[1], np.nan), np.full(unit.shape[1], np.nan)
    observed = unit.mean(axis=0)
    selected = unit[rng.integers(0, len(unit), size=(reps, len(unit)))].mean(axis=1)
    low, high = np.quantile(selected, [0.025, 0.975], axis=0)
    return observed, low, high


def slope_summary(
    records: list[dict[str, Any]], null: np.ndarray, scale: int, predictor: str,
    boot_reps: int, rng: np.random.Generator,
) -> list[dict[str, Any]]:
    indices = np.asarray([i for i, r in enumerate(records) if int(r["scale_m"]) == scale], int)
    birds = sorted({str(records[i]["individual"]) for i in indices})
    units: list[np.ndarray] = []; null_units: list[np.ndarray] = []; used_events = 0
    for bird in birds:
        use = np.asarray([i for i in indices if str(records[i]["individual"]) == bird], int)
        x = np.asarray([float(records[i][predictor]) for i in use], float)
        if len(use) < 3 or not np.isfinite(x).all() or np.var(x) <= 0:
            continue
        centered = x - x.mean(); denominator = float(centered @ centered)
        y = np.row_stack([records[i]["observed"] for i in use])
        units.append((centered[:, None] * y).sum(axis=0) / denominator)
        null_units.append((null[:, use, :] * centered[None, :, None]).sum(axis=1) / denominator)
        used_events += len(use)
    unit = np.row_stack(units) if units else np.empty((0, len(METRICS)))
    null_unit = np.stack(null_units, axis=1) if null_units else np.empty((len(null), 0, len(METRICS)))
    observed, low, high = bootstrap_units(unit, boot_reps, rng)
    null_population = null_unit.mean(axis=1) if len(units) else np.full((len(null), len(METRICS)), np.nan)
    rows = []
    for metric_index, metric in enumerate(METRICS):
        raw_p = float((1 + np.sum(null_population[:, metric_index] >= observed[metric_index])) / (len(null) + 1)) if len(units) else math.nan
        rows.append({
            "scale_m": scale, "predictor": predictor, "estimand": metric,
            "events": used_events, "individuals": len(units),
            "individual_equal_slope": observed[metric_index],
            "bootstrap_ci_low": low[metric_index], "bootstrap_ci_high": high[metric_index],
            "phase_null_mean": float(np.nanmean(null_population[:, metric_index])), "phase_p_raw": raw_p,
        })
    return rows


def group_summary(
    records: list[dict[str, Any]], null: np.ndarray, scale: int, group: str,
    boot_reps: int, rng: np.random.Generator,
) -> list[dict[str, Any]]:
    indices = np.asarray([
        i for i, r in enumerate(records) if int(r["scale_m"]) == scale and str(r["behavior_class"]) == group
    ], int)
    birds = sorted({str(records[i]["individual"]) for i in indices})
    units: list[np.ndarray] = []; null_units: list[np.ndarray] = []
    for bird in birds:
        use = np.asarray([i for i in indices if str(records[i]["individual"]) == bird], int)
        if not len(use):
            continue
        units.append(np.row_stack([records[i]["observed"] for i in use]).mean(axis=0))
        null_units.append(null[:, use, :].mean(axis=1))
    unit = np.row_stack(units) if units else np.empty((0, len(METRICS)))
    null_unit = np.stack(null_units, axis=1) if null_units else np.empty((len(null), 0, len(METRICS)))
    observed, low, high = bootstrap_units(unit, boot_reps, rng)
    null_population = null_unit.mean(axis=1) if len(units) else np.full((len(null), len(METRICS)), np.nan)
    rows = []
    for metric_index, metric in enumerate(METRICS):
        raw_p = float((1 + np.sum(null_population[:, metric_index] >= observed[metric_index])) / (len(null) + 1)) if len(units) else math.nan
        rows.append({
            "scale_m": scale, "behavior_group": group, "estimand": metric,
            "events": len(indices), "individuals": len(units),
            "individual_equal_mean": observed[metric_index],
            "bootstrap_ci_low": low[metric_index], "bootstrap_ci_high": high[metric_index],
            "phase_null_mean": float(np.nanmean(null_population[:, metric_index])), "phase_p_raw": raw_p,
        })
    return rows


def paired_summary(
    records: list[dict[str, Any]], null: np.ndarray, scale: int,
    boot_reps: int, rng: np.random.Generator,
) -> list[dict[str, Any]]:
    all_indices = [i for i, r in enumerate(records) if int(r["scale_m"]) == scale]
    birds = sorted({str(records[i]["individual"]) for i in all_indices})
    units: list[np.ndarray] = []; null_units: list[np.ndarray] = []; event_count = 0
    for bird in birds:
        forage = np.asarray([i for i in all_indices if str(records[i]["individual"]) == bird and records[i]["behavior_class"] == "forage_dominant"], int)
        rest = np.asarray([i for i in all_indices if str(records[i]["individual"]) == bird and records[i]["behavior_class"] == "rest_dominant"], int)
        if not len(forage) or not len(rest):
            continue
        units.append(
            np.row_stack([records[i]["observed"] for i in forage]).mean(axis=0)
            - np.row_stack([records[i]["observed"] for i in rest]).mean(axis=0)
        )
        null_units.append(null[:, forage, :].mean(axis=1) - null[:, rest, :].mean(axis=1))
        event_count += len(forage) + len(rest)
    unit = np.row_stack(units) if units else np.empty((0, len(METRICS)))
    null_unit = np.stack(null_units, axis=1) if null_units else np.empty((len(null), 0, len(METRICS)))
    observed, low, high = bootstrap_units(unit, boot_reps, rng)
    null_population = null_unit.mean(axis=1) if len(units) else np.full((len(null), len(METRICS)), np.nan)
    rows = []
    for metric_index, metric in enumerate(METRICS):
        raw_p = float((1 + np.sum(null_population[:, metric_index] >= observed[metric_index])) / (len(null) + 1)) if len(units) else math.nan
        rows.append({
            "scale_m": scale, "contrast": "forage_dominant_minus_rest_dominant", "estimand": metric,
            "events": event_count, "paired_individuals": len(units),
            "individual_equal_difference": observed[metric_index],
            "bootstrap_ci_low": low[metric_index], "bootstrap_ci_high": high[metric_index],
            "phase_null_mean": float(np.nanmean(null_population[:, metric_index])), "phase_p_raw": raw_p,
        })
    return rows


def add_holm(frame: pd.DataFrame, family_masks: dict[str, pd.Series]) -> pd.DataFrame:
    frame = frame.copy(); frame["holm_family"] = "descriptive"; frame["phase_p_holm"] = np.nan
    for family, mask in family_masks.items():
        index = frame.index[mask]
        if len(index):
            frame.loc[index, "holm_family"] = family
            frame.loc[index, "phase_p_holm"] = holm(frame.loc[index, "phase_p_raw"].tolist())
    return frame


def adjacent_pairs(scales: list[int]) -> list[tuple[int, int]]:
    ordered = sorted(scales)
    return [(a, b) for a, b in zip(ordered[:-1], ordered[1:])]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phase-reps", type=int, default=999)
    parser.add_argument("--bootstrap-reps", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 forbidden: affinity={affinity}")
    if args.phase_reps < 19 or args.bootstrap_reps < 100:
        raise RuntimeError("insufficient replicate counts")
    output = args.outdir; output.mkdir(parents=True, exist_ok=True)
    progress = Progress(output)
    rng = np.random.default_rng(args.seed)
    started = time.monotonic()

    cells, file_audit, chl_audit = audit_and_load_chl(progress)
    records, event_audit, behavior_input_audit = reconstruct_events(cells, progress)
    null, null_audit = phase_event_metrics(records, args.phase_reps, rng, progress)

    event_rows = []
    for record in records:
        row = {key: record[key] for key in (
            "dataset", "scale_m", "event_id", "run_key", "individual", "deployment",
            "origin_index_support_run", "endpoint_index_support_run", "future_specificity",
            "future_forage", "future_rest", "prior_specificity", "prior_forage", "prior_rest",
            "future_minus_prior_specificity", "behavior_class",
        )}
        row.update({name: float(value) for name, value in zip(METRICS, record["observed"])})
        row["n_native_cell_runs"] = len(record["values"])
        row["n_record_cell_runs"] = len(record["position"])
        event_rows.append(row)
    events = pd.DataFrame(event_rows)
    identity_error = float(np.max(np.abs(
        events[["E_high", "E_low", "E_union"]].to_numpy(float)
        - events[["R_high", "R_low", "R_union"]].to_numpy(float)
        - events[["L_high", "L_low", "L_union"]].to_numpy(float)
    )))
    if identity_error > 1e-12:
        raise RuntimeError("E=R+L identity failed")

    slope_rows = []
    predictors = (
        "future_specificity", "future_forage", "future_rest", "prior_specificity",
        "prior_forage", "prior_rest", "future_minus_prior_specificity",
    )
    for scale in SCALES:
        for predictor in predictors:
            slope_rows.extend(slope_summary(records, null, scale, predictor, args.bootstrap_reps, rng))
    slopes = pd.DataFrame(slope_rows)
    slopes = add_holm(slopes, {
        "primary_future_specificity_L_union_4_scales": slopes.predictor.eq("future_specificity") & slopes.estimand.eq("L_union"),
        "secondary_future_specificity_L_tails_8": slopes.predictor.eq("future_specificity") & slopes.estimand.isin(["L_high", "L_low"]),
        "secondary_future_specificity_E_R_union_8": slopes.predictor.eq("future_specificity") & slopes.estimand.isin(["E_union", "R_union"]),
    })

    group_rows = []
    for scale in SCALES:
        for group in ("forage_dominant", "rest_dominant"):
            group_rows.extend(group_summary(records, null, scale, group, args.bootstrap_reps, rng))
    groups = pd.DataFrame(group_rows)
    groups = add_holm(groups, {
        "primary_forage_dominant_L_union_4_scales": groups.behavior_group.eq("forage_dominant") & groups.estimand.eq("L_union"),
        "secondary_forage_dominant_L_tails_8": groups.behavior_group.eq("forage_dominant") & groups.estimand.isin(["L_high", "L_low"]),
    })

    paired_rows = []
    for scale in SCALES:
        paired_rows.extend(paired_summary(records, null, scale, args.bootstrap_reps, rng))
    paired = pd.DataFrame(paired_rows)
    paired = add_holm(paired, {
        "secondary_forage_minus_rest_L_union_4": paired.estimand.eq("L_union"),
        "secondary_forage_minus_rest_L_tails_8": paired.estimand.isin(["L_high", "L_low"]),
    })

    primary_slope = slopes.loc[
        slopes.predictor.eq("future_specificity") & slopes.estimand.eq("L_union")
    ].copy()
    primary_slope["coverage_pass"] = (primary_slope.events >= 30) & (primary_slope.individuals >= 6)
    primary_slope["pass"] = (
        primary_slope.coverage_pass & (primary_slope.individual_equal_slope > 0)
        & (primary_slope.bootstrap_ci_low > 0) & (primary_slope.phase_p_holm <= .05)
    )
    primary_group = groups.loc[
        groups.behavior_group.eq("forage_dominant") & groups.estimand.eq("L_union")
    ].copy()
    primary_group["coverage_pass"] = (primary_group.events >= 15) & (primary_group.individuals >= 5)
    primary_group["pass"] = (
        primary_group.coverage_pass & (primary_group.individual_equal_mean > 0)
        & (primary_group.bootstrap_ci_low > 0) & (primary_group.phase_p_holm <= .05)
    )
    combined_scales = sorted(set(primary_slope.loc[primary_slope["pass"], "scale_m"]) & set(primary_group.loc[primary_group["pass"], "scale_m"]))
    adjacent = [pair for pair in adjacent_pairs(list(SCALES)) if pair[0] in combined_scales and pair[1] in combined_scales]
    verdict = "DIRECT_FORAGING_STATE_LINK_SUPPORTED" if adjacent else "NO_DIRECT_FORAGING_STATE_LINK"

    events.to_csv(output / "event_metrics.csv.gz", index=False, compression="gzip")
    event_audit.to_csv(output / "event_eligibility_audit.csv", index=False)
    file_audit.to_csv(output / "chl_file_audit.csv", index=False)
    null_audit.to_csv(output / "phase_null_audit.csv", index=False)
    slopes.to_csv(output / "continuous_behavior_slopes.csv", index=False)
    groups.to_csv(output / "behavior_group_effects.csv", index=False)
    paired.to_csv(output / "paired_behavior_contrasts.csv", index=False)
    primary_slope.to_csv(output / "primary_slope_gate.csv", index=False)
    primary_group.to_csv(output / "primary_forage_group_gate.csv", index=False)

    sequence_path = output / "event_sequences.jsonl.gz"
    with gzip.open(sequence_path, "wt", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps({
                "scale_m": int(record["scale_m"]), "event_id": str(record["event_id"]),
                "run_key": str(record["run_key"]), "individual": str(record["individual"]),
                "values": finite(record["values"]), "record_runs": finite(record["position"]),
            }, ensure_ascii=False) + "\n")

    summary = {
        "status": verdict, "formal_seed": args.seed, "phase_reps": args.phase_reps,
        "bootstrap_reps": args.bootstrap_reps, "cpu_affinity": affinity,
        "runtime_s": time.monotonic() - started, "eligible_events": len(events),
        "eligible_individuals": int(events.individual.nunique()),
        "events_by_scale": {str(k): int(v) for k, v in events.groupby("scale_m").size().items()},
        "individuals_by_scale": {str(k): int(v) for k, v in events.groupby("scale_m").individual.nunique().items()},
        "maximum_E_minus_R_minus_L_error": identity_error,
        "primary_slope_pass_scales": primary_slope.loc[primary_slope["pass"], "scale_m"].astype(int).tolist(),
        "primary_forage_group_pass_scales": primary_group.loc[primary_group["pass"], "scale_m"].astype(int).tolist(),
        "combined_pass_scales": [int(v) for v in combined_scales],
        "adjacent_combined_pass_pairs": [[int(a), int(b)] for a, b in adjacent],
        "chl_audit": chl_audit, "behavior_input_audit": behavior_input_audit,
        "hashes": {
            "prereg": sha256(PREREG), "catalog": sha256(CATALOG), "manifest": sha256(MANIFEST),
            "frozen_endpoints": sha256(ENDPOINTS), "public_behavior": sha256(short.DATA),
            "script": sha256(Path(__file__)), "event_metrics": sha256(output / "event_metrics.csv.gz"),
            "event_sequences": sha256(sequence_path),
        },
    }
    (output / "summary.json").write_text(json.dumps(finite(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = [
        "# D107 short-tailed shearwater last-passage CHL × behavior result", "",
        f"- Verdict: `{verdict}`", f"- Eligible events: {len(events)}; individuals: {events.individual.nunique()}",
        f"- Primary slope pass scales: {summary['primary_slope_pass_scales']}",
        f"- Primary forage-group pass scales: {summary['primary_forage_group_pass_scales']}",
        f"- Adjacent joint-pass pairs: {summary['adjacent_combined_pass_pairs']}", "",
        "## Primary continuous interaction", "", primary_slope.to_markdown(index=False), "",
        "## Primary forage-dominant within-group footprint", "", primary_group.to_markdown(index=False), "",
        "## Interpretation boundary", "",
        "A positive gate supports a same-event link between the CHL last-passage footprint and foraging-versus-rest state. It does not establish direct CHL perception, a sensory channel, successful prey capture, or causal triggering.", "",
    ]
    (output / "report.md").write_text("\n".join(report), encoding="utf-8")
    progress.emit("complete", 1, 1, verdict=verdict)
    print(json.dumps(finite(summary), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
