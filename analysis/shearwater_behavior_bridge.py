#!/usr/bin/env python3
"""Five-minute-support RD bridge to short-tailed shearwater accelerometry states."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import CRS, Transformer


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA = ROOT / "external_inputs/shearwater_behavior/GPS_tracking_and_behavioural_states.csv"
DEFAULT_OUTPUT = ROOT / "results/shearwater_behavior_bridge"
PHASES_S = (0, 60, 120, 180, 240)
DELTAS_M = (500.0, 1000.0, 2000.0, 5000.0, 10000.0)
WINDOWS = ("FUTURE_5MIN", "PRIOR_5MIN")
METRICS = ("specificity", "forage", "rest", "flight", "takeoff")
BEHAVIOR_MAP = {
    "forage": {"diving", "sforaging"},
    "rest": {"resting"},
    "flight": {"flapping", "gliding"},
    "takeoff": {"takeoff"},
}

from radial_drawdown.core import detect as detect_sampled_fix  # noqa: E402
from paper2_core import extract_continuous_drawdown_events_nd  # noqa: E402


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


class Progress:
    def __init__(self, path: Path, total: int):
        self.path = path
        self.total = total
        self.started = time.time()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("", encoding="utf-8")

    def emit(self, stage: str, done: int, **extra: object) -> None:
        elapsed = max(time.time() - self.started, 1e-9)
        rate = done / elapsed if done else 0.0
        eta = (self.total - done) / rate if rate and done <= self.total else None
        row = {
            "stage": stage, "completed": done, "total": self.total,
            "percent": round(100.0 * done / self.total, 1),
            "elapsed_s": round(elapsed, 3), "throughput_stages_per_s": round(rate, 4),
            "eta_s": None if eta is None else round(eta, 3), **extra,
        }
        line = json.dumps(row, ensure_ascii=False)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
        print(line, flush=True)


def individual_from_deployment(value: str) -> str:
    return value.split("-(", 1)[0]


def rolling_fraction(state: np.ndarray, labels: set[str], direction: str) -> np.ndarray:
    z = np.isin(state, list(labels)).astype(np.int64)
    cs = np.r_[0, np.cumsum(z)]
    out = np.full(len(z), np.nan)
    if direction == "future":
        idx = np.arange(0, max(0, len(z) - 299))
        out[idx] = (cs[idx + 300] - cs[idx]) / 300.0
    elif direction == "prior":
        idx = np.arange(299, len(z))
        out[idx] = (cs[idx + 1] - cs[idx - 299]) / 300.0
    else:
        raise KeyError(direction)
    return out


def load_highres(max_deployments: int) -> tuple[list[pd.DataFrame], pd.DataFrame, dict[str, object]]:
    d = pd.read_csv(DATA, dtype={"ID": str}, low_memory=False)
    d["timestamp"] = pd.to_datetime(d["Date"], format="%d/%m/%y %H:%M:%S", errors="coerce")
    d["lat"] = pd.to_numeric(d["Latitude"], errors="coerce")
    d["lon"] = pd.to_numeric(d["Longitude"], errors="coerce")
    d["deployment"] = d["ID"].astype(str)
    d["individual"] = d["deployment"].map(individual_from_deployment)
    d["behavior"] = d["Behaviour"].astype(str)
    d["_source_row"] = np.arange(len(d), dtype=np.int64)
    input_rows = len(d)
    deployments = sorted(d["deployment"].dropna().unique())
    if max_deployments > 0:
        deployments = deployments[:max_deployments]
        d = d[d["deployment"].isin(deployments)].copy()
    d = d[
        d["timestamp"].notna() & d["lat"].between(-90, 90) & d["lon"].between(-180, 180)
        & d["behavior"].isin(set().union(*BEHAVIOR_MAP.values()))
    ].copy()
    d = d.sort_values(["deployment", "timestamp", "_source_row"], kind="stable")
    before_duplicate = len(d)
    d = d.drop_duplicates(["deployment", "timestamp"], keep="first").copy()
    duplicates_removed = before_duplicate - len(d)

    center_lat = float(d["lat"].median())
    center_lon = float(d["lon"].median())
    crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={center_lat:.12f} +lon_0={center_lon:.12f} +datum=WGS84 +units=m +no_defs"
    )
    transformer = Transformer.from_crs(4326, crs, always_xy=True)
    d["x_m"], d["y_m"] = transformer.transform(d["lon"].to_numpy(float), d["lat"].to_numpy(float))

    high_runs: list[pd.DataFrame] = []
    deployment_rows: list[dict[str, object]] = []
    total_non1s = 0
    total_impossible = 0
    for deployment, g0 in d.groupby("deployment", sort=True):
        g = g0.sort_values(["timestamp", "_source_row"], kind="stable").reset_index(drop=True)
        dt = g["timestamp"].diff().dt.total_seconds().to_numpy(float)
        distance = np.hypot(g["x_m"].diff().to_numpy(float), g["y_m"].diff().to_numpy(float))
        speed = np.divide(distance, dt, out=np.full(len(g), np.nan), where=dt > 0)
        source_adjacent = g["_source_row"].diff().eq(1).to_numpy()
        edge_ok = (dt == 1.0) & source_adjacent & (~np.isfinite(speed) | (speed <= 100.0))
        new_run = ~edge_ok
        new_run[0] = True
        g["local_run"] = np.cumsum(new_run)
        runs_this = 0
        points_this = 0
        for local_run, r0 in g.groupby("local_run", sort=False):
            r = r0.reset_index(drop=True).copy()
            if len(r) < 900:  # at least three 5-min supports for every useful phase
                continue
            state = r["behavior"].to_numpy(str)
            for name, labels in BEHAVIOR_MAP.items():
                r[f"FUTURE_5MIN__{name}"] = rolling_fraction(state, labels, "future")
                r[f"PRIOR_5MIN__{name}"] = rolling_fraction(state, labels, "prior")
            r["FUTURE_5MIN__specificity"] = r["FUTURE_5MIN__forage"] - r["FUTURE_5MIN__rest"]
            r["PRIOR_5MIN__specificity"] = r["PRIOR_5MIN__forage"] - r["PRIOR_5MIN__rest"]
            r["high_run_key"] = f"{deployment}:H{int(local_run):03d}"
            high_runs.append(r)
            runs_this += 1
            points_this += len(r)
        non1s = int(np.sum(np.isfinite(dt) & (dt != 1.0)))
        impossible = int(np.sum(np.isfinite(speed) & (speed > 100.0)))
        total_non1s += non1s
        total_impossible += impossible
        deployment_rows.append({
            "deployment": deployment, "individual": individual_from_deployment(deployment),
            "valid_unique_rows": len(g), "eligible_highres_runs_ge900": runs_this,
            "eligible_highres_points": points_this, "non_1s_edges": non1s,
            "impossible_speed_edges": impossible,
            "behavior_forage_fraction": float(g["behavior"].isin(BEHAVIOR_MAP["forage"]).mean()),
            "behavior_rest_fraction": float(g["behavior"].isin(BEHAVIOR_MAP["rest"]).mean()),
        })
    audit = {
        "input_rows_all_deployments": input_rows,
        "selected_valid_unique_rows": len(d),
        "deployments": int(d["deployment"].nunique()),
        "biological_individuals": int(d["individual"].nunique()),
        "duplicates_removed": duplicates_removed,
        "non_1s_edges": total_non1s,
        "impossible_speed_edges": total_impossible,
        "eligible_highres_runs_ge900": len(high_runs),
        "eligible_highres_points": int(sum(len(r) for r in high_runs)),
        "hardware_gps_interval_s": 300,
        "public_behavior_interval_s": 1,
        "spatial_interpretation": "per-second positions are synchronization/interpolated support, not independent 1-Hz GPS fixes",
        "projection": "WGS84_AEQD",
        "projection_center_lat": center_lat,
        "projection_center_lon": center_lon,
    }
    return high_runs, pd.DataFrame(deployment_rows), audit


def build_support_runs(high_runs: list[pd.DataFrame]) -> dict[int, list[pd.DataFrame]]:
    output: dict[int, list[pd.DataFrame]] = {phase: [] for phase in PHASES_S}
    for high in high_runs:
        for phase in PHASES_S:
            idx = np.arange(phase, len(high), 300, dtype=int)
            support = high.iloc[idx].copy().reset_index(drop=True)
            if len(support) < 3:
                continue
            support["support_source_index_highrun"] = idx
            support["elapsed_s"] = (
                support["timestamp"] - support["timestamp"].iloc[0]
            ).dt.total_seconds()
            support["phase_s"] = phase
            support["run_key"] = f"{high['high_run_key'].iloc[0]}:P{phase:03d}"
            output[phase].append(support)
    return output


def endpoint_masks(support_by_phase: dict[int, list[pd.DataFrame]]) -> tuple[dict[tuple[int, str, float, str], np.ndarray], pd.DataFrame]:
    masks: dict[tuple[int, str, float, str], np.ndarray] = {}
    rows: list[dict[str, object]] = []
    for phase, runs in support_by_phase.items():
        for run in runs:
            key = str(run["run_key"].iloc[0])
            xy = run[["x_m", "y_m"]].to_numpy(float)
            seconds = run["elapsed_s"].to_numpy(float)
            for delta in DELTAS_M:
                pmask = np.zeros(len(run), dtype=bool)
                for record in detect_sampled_fix(seconds, xy, delta).records():
                    if record["censored"]:
                        continue
                    idx = int(record["endpoint_index"])
                    pmask[idx] = True
                    rows.append({
                        "phase_s": phase, "implementation": "RD_P_SAMPLED_FIX", "delta_m": int(delta),
                        "individual": str(run["individual"].iloc[0]), "deployment": str(run["deployment"].iloc[0]),
                        "run_key": key, "endpoint_index_support_run": idx,
                        "endpoint_timestamp": run.loc[idx, "timestamp"], "endpoint_source_row": int(run.loc[idx, "_source_row"]),
                    })
                masks[(phase, "RD_P_SAMPLED_FIX", delta, key)] = pmask

                lmask = np.zeros(len(run), dtype=bool)
                for event in extract_continuous_drawdown_events_nd(xy, delta):
                    if event.is_terminal:
                        continue
                    idx = int(event.endpoint_idx)
                    lmask[idx] = True
                    rows.append({
                        "phase_s": phase, "implementation": "RD_L_CONTINUOUS_LINE", "delta_m": int(delta),
                        "individual": str(run["individual"].iloc[0]), "deployment": str(run["deployment"].iloc[0]),
                        "run_key": key, "endpoint_index_support_run": idx,
                        "endpoint_timestamp": run.loc[idx, "timestamp"], "endpoint_source_row": int(run.loc[idx, "_source_row"]),
                    })
                masks[(phase, "RD_L_CONTINUOUS_LINE", delta, key)] = lmask
    return masks, pd.DataFrame(rows)


def circular_weighted_overlap(selector: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.fft.ifft(np.conj(np.fft.fft(selector.astype(float))) * np.fft.fft(values.astype(float))).real


def bootstrap_ci(values: np.ndarray, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    if len(values) == 0:
        return math.nan, math.nan
    boot = rng.choice(values, size=(n_boot, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(boot, .025)), float(np.quantile(boot, .975))


def association_test(
    runs: list[pd.DataFrame], masks: dict[tuple[int, str, float, str], np.ndarray], phase: int,
    implementation: str, delta: float, window: str, metric: str, n_perm: int, n_boot: int, seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    individuals = sorted({str(r["individual"].iloc[0]) for r in runs})
    endpoint_n = {bird: 0 for bird in individuals}
    observed_sum = {bird: 0.0 for bird in individuals}
    null_sum = {bird: np.zeros(n_perm, dtype=float) for bird in individuals}
    used_runs = 0
    for run in runs:
        key = str(run["run_key"].iloc[0])
        bird = str(run["individual"].iloc[0])
        values = run[f"{window}__{metric}"].to_numpy(float)
        valid = np.isfinite(values)
        if int(valid.sum()) < 3:
            continue
        endpoint = masks[(phase, implementation, delta, key)][valid]
        values = values[valid]
        n_endpoint = int(endpoint.sum())
        if n_endpoint == 0:
            continue
        used_runs += 1
        endpoint_n[bird] += n_endpoint
        observed_sum[bird] += float(np.sum(values[endpoint]))
        sums_by_shift = circular_weighted_overlap(endpoint, values)
        shifts = rng.integers(0, len(values), size=n_perm)
        null_sum[bird] += sums_by_shift[shifts]

    eligible = [bird for bird in individuals if endpoint_n[bird] > 0]
    if eligible:
        observed_by = np.asarray([observed_sum[b] / endpoint_n[b] for b in eligible], dtype=float)
        null_by = np.vstack([null_sum[b] / endpoint_n[b] for b in eligible])
        expected_by = null_by.mean(axis=1)
        excess_by = observed_by - expected_by
        observed_equal = float(observed_by.mean())
        null_equal = null_by.mean(axis=0)
        excess = float(excess_by.mean())
        p = float((1 + np.sum(null_equal >= observed_equal - 1e-15)) / (n_perm + 1))
        ci_low, ci_high = bootstrap_ci(excess_by, n_boot, np.random.default_rng(seed + 1000003))
    else:
        observed_by = expected_by = excess_by = np.asarray([], dtype=float)
        observed_equal = excess = p = ci_low = ci_high = math.nan
    by_individual = pd.DataFrame({
        "individual": eligible,
        "endpoints": [endpoint_n[b] for b in eligible],
        "observed_endpoint_mean": observed_by,
        "circular_null_expected_mean": expected_by,
        "endpoint_excess": excess_by,
    })
    result = {
        "phase_s": phase, "implementation": implementation, "delta_m": int(delta),
        "window": window, "metric": metric, "n_individuals": len(eligible),
        "n_endpoints": int(sum(endpoint_n[b] for b in eligible)), "n_runs_with_endpoint": used_runs,
        "individual_equal_observed_endpoint_mean": observed_equal,
        "individual_equal_circular_null_mean": float(observed_equal - excess) if np.isfinite(excess) else math.nan,
        "individual_equal_endpoint_excess": excess,
        "individual_bootstrap_ci_low": ci_low, "individual_bootstrap_ci_high": ci_high,
        "circular_shift_p_one_sided": p,
        "individuals_positive": int(np.sum(excess_by > 0)),
        "individuals_negative": int(np.sum(excess_by < 0)),
        "individuals_zero": int(np.sum(excess_by == 0)),
    }
    return result, by_individual


def holm(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    order = np.argsort(np.where(np.isfinite(p), p, 1.0))
    adjusted = np.full(len(p), np.nan)
    running = 0.0
    for rank, idx in enumerate(order):
        raw = p[idx] if np.isfinite(p[idx]) else 1.0
        running = max(running, min(1.0, (len(p) - rank) * raw))
        adjusted[idx] = running
    return adjusted.tolist()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-perm", type=int, default=10000)
    parser.add_argument("--n-boot", type=int, default=10000)
    parser.add_argument("--max-deployments", type=int, default=0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.output / "progress.jsonl", 7)
    started = time.time()

    high_runs, deployment_audit, input_audit = load_highres(args.max_deployments)
    deployment_audit.to_csv(args.output / "deployment_and_behavior_audit.csv", index=False)
    (args.output / "input_audit.json").write_text(json.dumps(input_audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("input_and_resolution_audit", 1, high_runs=len(high_runs), highres_points=sum(len(r) for r in high_runs), deployments=input_audit["deployments"], individuals=input_audit["biological_individuals"])

    support_by_phase = build_support_runs(high_runs)
    support_audit = pd.DataFrame([
        {"phase_s": phase, "support_runs": len(runs), "support_points": int(sum(len(r) for r in runs)), "individuals": len({str(r['individual'].iloc[0]) for r in runs})}
        for phase, runs in support_by_phase.items()
    ])
    support_audit.to_csv(args.output / "five_minute_support_audit.csv", index=False)
    progress.emit("five_minute_support", 2, phase0_points=int(support_audit.loc[support_audit["phase_s"].eq(0), "support_points"].iloc[0]))

    masks, endpoints = endpoint_masks(support_by_phase)
    endpoints.to_csv(args.output / "behavior_blind_five_minute_rd_endpoints.csv", index=False)
    progress.emit("behavior_blind_rd", 3, endpoints=len(endpoints))

    results: list[dict[str, object]] = []
    individual_frames: list[pd.DataFrame] = []
    test_index = 0
    for phase in PHASES_S:
        for implementation in ("RD_P_SAMPLED_FIX", "RD_L_CONTINUOUS_LINE"):
            for delta in DELTAS_M:
                for window in WINDOWS:
                    for metric in METRICS:
                        test_index += 1
                        result, by_individual = association_test(
                            support_by_phase[phase], masks, phase, implementation, delta, window, metric,
                            args.n_perm, args.n_boot, 2026082300 + test_index * 1009,
                        )
                        results.append(result)
                        if not by_individual.empty:
                            by_individual.insert(0, "metric", metric)
                            by_individual.insert(0, "window", window)
                            by_individual.insert(0, "delta_m", int(delta))
                            by_individual.insert(0, "implementation", implementation)
                            by_individual.insert(0, "phase_s", phase)
                            individual_frames.append(by_individual)
    association = pd.DataFrame(results)
    association.to_csv(args.output / "all_rd_behavior_associations_preholm.csv", index=False)
    pd.concat(individual_frames, ignore_index=True).to_csv(args.output / "all_rd_behavior_associations_by_individual.csv", index=False)
    progress.emit("association_tests", 4, tests=len(association))

    phase0_future = association[
        association["phase_s"].eq(0) & association["implementation"].eq("RD_P_SAMPLED_FIX")
        & association["window"].eq("FUTURE_5MIN") & association["metric"].isin(["specificity", "forage"])
    ].copy()
    phase0_future["holm_p_within_metric_five_scales"] = np.nan
    for metric in ("specificity", "forage"):
        idx = phase0_future.index[phase0_future["metric"].eq(metric)]
        phase0_future.loc[idx, "holm_p_within_metric_five_scales"] = holm(phase0_future.loc[idx, "circular_shift_p_one_sided"].tolist())
    association = association.merge(
        phase0_future[["phase_s", "implementation", "delta_m", "window", "metric", "holm_p_within_metric_five_scales"]],
        on=["phase_s", "implementation", "delta_m", "window", "metric"], how="left",
    )
    association.to_csv(args.output / "all_rd_behavior_associations.csv", index=False)

    primary = association[
        association["phase_s"].eq(0) & association["implementation"].eq("RD_P_SAMPLED_FIX")
        & association["window"].eq("FUTURE_5MIN")
    ].copy()
    primary_pivot = primary.pivot(index="delta_m", columns="metric", values="individual_equal_endpoint_excess")
    scale_rows: list[dict[str, object]] = []
    for delta in DELTAS_M:
        spec = primary[(primary["delta_m"].eq(delta)) & primary["metric"].eq("specificity")].iloc[0]
        forage = primary[(primary["delta_m"].eq(delta)) & primary["metric"].eq("forage")].iloc[0]
        rest = primary[(primary["delta_m"].eq(delta)) & primary["metric"].eq("rest")].iloc[0]
        phase_signs = []
        for phase in PHASES_S:
            z = association[
                association["phase_s"].eq(phase) & association["implementation"].eq("RD_P_SAMPLED_FIX")
                & association["delta_m"].eq(delta) & association["window"].eq("FUTURE_5MIN")
                & association["metric"].eq("specificity")
            ].iloc[0]
            phase_signs.append(float(z["individual_equal_endpoint_excess"]) > 0)
        rd_l = association[
            association["phase_s"].eq(0) & association["implementation"].eq("RD_L_CONTINUOUS_LINE")
            & association["delta_m"].eq(delta) & association["window"].eq("FUTURE_5MIN")
            & association["metric"].eq("specificity")
        ].iloc[0]
        gate = bool(
            spec["n_individuals"] >= 8 and spec["n_endpoints"] >= 30
            and spec["individual_equal_endpoint_excess"] > 0 and spec["individual_bootstrap_ci_low"] > 0
            and spec["holm_p_within_metric_five_scales"] < .05
            and forage["individual_equal_endpoint_excess"] > 0 and forage["individual_bootstrap_ci_low"] > 0
            and forage["individual_equal_endpoint_excess"] > rest["individual_equal_endpoint_excess"]
            and sum(phase_signs) >= 4 and rd_l["individual_equal_endpoint_excess"] >= 0
        )
        scale_rows.append({
            "delta_m": int(delta), "n_individuals": int(spec["n_individuals"]), "n_endpoints": int(spec["n_endpoints"]),
            "specificity_excess": spec["individual_equal_endpoint_excess"],
            "specificity_ci_low": spec["individual_bootstrap_ci_low"], "specificity_ci_high": spec["individual_bootstrap_ci_high"],
            "specificity_holm_p": spec["holm_p_within_metric_five_scales"],
            "forage_excess": forage["individual_equal_endpoint_excess"],
            "forage_ci_low": forage["individual_bootstrap_ci_low"], "forage_ci_high": forage["individual_bootstrap_ci_high"],
            "rest_excess": rest["individual_equal_endpoint_excess"],
            "positive_phases_of_5": int(sum(phase_signs)), "rd_l_specificity_excess": rd_l["individual_equal_endpoint_excess"],
            "scale_gate": gate,
        })
    scale_family = pd.DataFrame(scale_rows)
    scale_family.to_csv(args.output / "primary_scale_family.csv", index=False)
    gates = scale_family["scale_gate"].to_numpy(bool)
    adjacent_pass = any(gates[i] and gates[i + 1] for i in range(len(gates) - 1))
    primary_pass = bool(adjacent_pass)
    verdict = "FORAGING_SPECIFIC_RENEWAL_BRIDGE" if primary_pass else "NO_FORAGING_SPECIFIC_RENEWAL_BRIDGE_PASS"
    progress.emit("primary_verdict", 5, verdict=verdict, passing_scales=scale_family.loc[scale_family["scale_gate"], "delta_m"].astype(int).tolist())

    coverage = endpoints.groupby(["phase_s", "implementation", "delta_m"], as_index=False).agg(endpoints=("endpoint_index_support_run", "size"), individuals=("individual", "nunique"), deployments=("deployment", "nunique"))
    coverage.to_csv(args.output / "endpoint_coverage_summary.csv", index=False)
    progress.emit("coverage", 6, rows=len(coverage))

    prior_primary = association[
        association["phase_s"].eq(0) & association["implementation"].eq("RD_P_SAMPLED_FIX")
        & association["window"].eq("PRIOR_5MIN") & association["metric"].isin(["specificity", "forage", "rest"])
    ][["delta_m", "metric", "n_individuals", "n_endpoints", "individual_equal_endpoint_excess", "individual_bootstrap_ci_low", "individual_bootstrap_ci_high", "circular_shift_p_one_sided"]]
    report = f"""# 短尾鹱foraging-vs-rest行为特异性桥结果

正式裁决：`{verdict}`。

输入为{input_audit['deployments']}个deployment/{input_audit['biological_individuals']}个生物individual；原始硬件GPS时间支持为5 min，逐秒行为来自加速度分类。phase 0有{int(support_audit.loc[support_audit['phase_s'].eq(0), 'support_points'].iloc[0])}个规则5-min空间支持点。

## 主尺度族（phase 0、RD-P、未来5 min）

{scale_family.to_markdown(index=False)}

## endpoint前5 min方向诊断

{prior_primary.to_markdown(index=False)}

## 五相位/实现覆盖

{coverage.to_markdown(index=False)}

## 解释边界

- 逐秒坐标不是独立1-Hz GPS；本分析只使用规则5-min空间支持并报告五相位。
- `diving`没有压力传感器或捕食成功真值；foraging为加速度推断行为状态。
- 阳性最多支持renewal与foraging相对rest的状态特异关联，不识别cue、感觉通道、CHL直接感知或Lévy生成。
"""
    (args.output / "formal_report_CN.md").write_text(report, encoding="utf-8")
    formal = {
        "verdict": verdict, "primary_pass": primary_pass,
        "passing_scales_m": scale_family.loc[scale_family["scale_gate"], "delta_m"].astype(int).tolist(),
        "n_deployments": input_audit["deployments"], "n_individuals": input_audit["biological_individuals"],
        "phase0_support_points": int(support_audit.loc[support_audit["phase_s"].eq(0), "support_points"].iloc[0]),
        "n_permutations": args.n_perm, "n_bootstrap": args.n_boot, "elapsed_s": time.time() - started,
        "input_sha256": digest(DATA), "paper1_core_sha256": digest(PAPER1_CORE),
        "resolution_correction": "5-min GPS hardware support; per-second coordinates not treated as independent GPS",
    }
    (args.output / "formal_result.json").write_text(json.dumps(formal, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    progress.emit("complete", 7, elapsed_s=round(time.time() - started, 3))


if __name__ == "__main__":
    main()
