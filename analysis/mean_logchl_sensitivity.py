#!/usr/bin/env python3
"""D130: one frozen mean(log CHL) sensitivity for the D127 joint model.

This analysis reuses the exact D127 events and last-record metrics. It does not
rerun radial drawdown, environmental annotation, event reconstruction or phase
randomisation. The only model change is replacement of endpoint-excluded
median(log CHL) by endpoint-excluded mean(log CHL).
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
sys.path.insert(0, str(ROOT / "analysis"))
import laysan_same_event_integration as d127  # noqa: E402

PREREG = ROOT / "metadata/mean_logchl_sensitivity_specification_cn.md"
D43 = ROOT / "external_inputs/laysan_same_event/absolute_chl_step_event_table.csv.gz"
D104_A = ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz"
D127_OUT = ROOT / "results/laysan_same_event"
DEFAULT_OUT = ROOT / "results/mean_logchl_sensitivity"
SCALES = (500, 1000, 2000)
SEED = 1300827


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
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


def prepare_model(frame: pd.DataFrame) -> pd.DataFrame:
    columns = ["step_length_km", "interior_mean_logchl", "L_low", "L_high", "L_union", "cluster"]
    d = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=columns).copy()
    d = d[d.step_length_km > 0].copy()
    d["log_length"] = np.log(d.step_length_km)
    d["z_log_length"] = d127.within_cluster_z(d, "log_length")
    d["z_absolute"] = d127.within_cluster_z(d, "interior_mean_logchl")
    d["L_low_dm"] = d.L_low - d.groupby("cluster").L_low.transform("mean")
    d["L_high_dm"] = d.L_high - d.groupby("cluster").L_high.transform("mean")
    for source, target in (("interior_mean_logchl", "abs_tertile"), ("log_length", "length_tertile")):
        rank = d.groupby("cluster")[source].rank(method="average", pct=True)
        d[target] = np.minimum((rank * 3).astype(int), 2)
    return d


def bootstrap_summary(frame: pd.DataFrame, reps: int, rng: np.random.Generator) -> tuple[dict[str, Any], pd.DataFrame]:
    d = prepare_model(frame)
    clusters = d.cluster.drop_duplicates().to_numpy(object)
    by = {cluster: d[d.cluster == cluster] for cluster in clusters}
    beta = d127.fit_fe(d)
    animal_grid = d.groupby(["cluster", "abs_tertile", "length_tertile"], as_index=False).L_union.mean()
    animal_stat = animal_grid.groupby("cluster").L_union.mean()
    estimate_l = float(animal_stat.mean())
    grid = animal_grid.groupby(["abs_tertile", "length_tertile"], as_index=False).L_union.mean()

    xtx = np.empty((len(clusters), 3, 3), float)
    xty = np.empty((len(clusters), 3), float)
    l_values = np.empty(len(clusters), float)
    for j, cluster in enumerate(clusters):
        part = by[cluster]
        x = part[["z_absolute", "L_low_dm", "L_high_dm"]].to_numpy(float)
        y = part["z_log_length"].to_numpy(float)
        xtx[j] = x.T @ x
        xty[j] = x.T @ y
        l_values[j] = float(animal_stat.loc[cluster])

    choices = rng.integers(0, len(clusters), size=(reps, len(clusters)))
    counts = np.zeros((reps, len(clusters)), dtype=np.int16)
    np.add.at(counts, (np.repeat(np.arange(reps), len(clusters)), choices.ravel()), 1)
    boot_xtx = np.einsum("rc,cij->rij", counts, xtx, optimize=True)
    boot_xty = np.einsum("rc,cj->rj", counts, xty, optimize=True)
    beta_boot = np.einsum("rij,rj->ri", np.linalg.pinv(boot_xtx), boot_xty, optimize=True)
    l_boot = counts @ l_values / len(clusters)
    beta_low, beta_high = np.quantile(beta_boot, [0.025, 0.975], axis=0)
    l_low, l_high = np.quantile(l_boot, [0.025, 0.975])

    summary = {
        "events": len(d),
        "animals": len(clusters),
        "beta_absolute": beta[0],
        "beta_absolute_ci_low": beta_low[0],
        "beta_absolute_ci_high": beta_high[0],
        "beta_L_low": beta[1],
        "beta_L_low_ci_low": beta_low[1],
        "beta_L_low_ci_high": beta_high[1],
        "beta_L_high": beta[2],
        "beta_L_high_ci_low": beta_low[2],
        "beta_L_high_ci_high": beta_high[2],
        "standardized_L_union": estimate_l,
        "standardized_L_union_ci_low": l_low,
        "standardized_L_union_ci_high": l_high,
        "positive_grids": int((grid.L_union > 0).sum()),
        "estimable_grids": len(grid),
        "bootstrap_reps": reps,
    }
    return summary, grid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bootstrap-reps", type=int, default=20000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 forbidden: {affinity}")
    started = time.monotonic()
    print(f"PROGRESS start 0/4 (0.0%) affinity={affinity}", flush=True)

    d43_columns = [
        "dataset", "scale_m", "event_id", "cluster_id", "step_length_km",
        "interior_median_logchl", "interior_mean_logchl",
    ]
    d43 = pd.read_csv(D43, usecols=d43_columns, low_memory=False)
    d104a = pd.read_csv(D104_A, low_memory=False)
    goto = d104a[d104a.dataset.eq("goto")].merge(
        d43[d43.dataset.eq("goto")], on=["dataset", "scale_m", "event_id"], validate="one_to_one"
    )
    goto["cluster"] = goto["cluster"].astype(str)

    laysan_metrics = pd.read_csv(D127_OUT / "laysan_last_passage_event_metrics.csv.gz", low_memory=False)
    laysan_background = d43[d43.dataset.eq("usgs_laysan_albatross") & d43.scale_m.isin(SCALES)].copy()
    laysan = laysan_metrics.merge(
        laysan_background[
            ["dataset", "scale_m", "event_id", "step_length_km", "interior_median_logchl", "interior_mean_logchl"]
        ],
        on=["dataset", "scale_m", "event_id"],
        validate="one_to_one",
    )
    laysan["cluster"] = laysan["cluster"].astype(str)

    rng = np.random.default_rng(SEED)
    rows: list[dict[str, Any]] = []
    grids: list[pd.DataFrame] = []
    systems = [("goto", 100, goto)] + [
        ("usgs_laysan_albatross", scale, laysan[laysan.scale_m.eq(scale)]) for scale in SCALES
    ]
    for done, (dataset, scale, frame) in enumerate(systems, 1):
        summary, grid = bootstrap_summary(frame, args.bootstrap_reps, rng)
        summary.update({"dataset": dataset, "scale_m": scale})
        rows.append(summary)
        grid.insert(0, "scale_m", scale)
        grid.insert(0, "dataset", dataset)
        grids.append(grid)
        elapsed = time.monotonic() - started
        print(f"PROGRESS joint_model {done}/4 ({25 * done:.1f}%) elapsed={elapsed:.1f}s", flush=True)

    model = pd.DataFrame(rows)
    model["conditional_L_pass"] = model.standardized_L_union_ci_low.gt(0) & model.positive_grids.ge(6)
    model["absolute_pass"] = model.beta_absolute_ci_high.lt(0)
    model["L_low_positive"] = model.beta_L_low.gt(0)
    model["L_low_strict"] = model.beta_L_low_ci_low.gt(0)
    laysan_pass = sorted(
        model.loc[model.dataset.eq("usgs_laysan_albatross") & model.conditional_L_pass & model.absolute_pass, "scale_m"]
        .astype(int)
        .tolist()
    )
    adjacent = any(left in laysan_pass and right in laysan_pass for left, right in zip(SCALES[:-1], SCALES[1:]))
    goto_row = model[model.dataset.eq("goto")].iloc[0]
    l_low_strict_cross = bool(
        goto_row.L_low_strict and model.loc[model.dataset.eq("usgs_laysan_albatross"), "L_low_strict"].any()
    )
    all_passing_low_positive = bool(
        model.loc[model.conditional_L_pass & model.absolute_pass, "L_low_positive"].all()
    )
    gate = bool(
        goto_row.conditional_L_pass
        and goto_row.absolute_pass
        and adjacent
        and l_low_strict_cross
        and all_passing_low_positive
    )

    grid_table = pd.concat(grids, ignore_index=True)
    model.to_csv(args.output / "mean_logchl_joint_model_summary.csv", index=False)
    grid_table.to_csv(args.output / "mean_logchl_conditional_L_3x3_grid.csv", index=False)
    comparison = model.merge(
        pd.read_csv(D127_OUT / "joint_model_summary.csv"),
        on=["dataset", "scale_m"],
        suffixes=("_mean", "_median"),
        validate="one_to_one",
    )
    comparison.to_csv(args.output / "mean_vs_median_joint_model_comparison.csv", index=False)
    payload = {
        "status": "FORMAL_COMPLETE",
        "verdict": "MEAN_LOGCHL_SENSITIVITY_PASSED" if gate else "MEAN_LOGCHL_SENSITIVITY_FAILED",
        "gate_pass": gate,
        "laysan_scales_passing_absolute_and_conditional": laysan_pass,
        "laysan_adjacent_scale_pass": adjacent,
        "cross_system_L_low_strict": l_low_strict_cross,
        "all_passing_systems_L_low_positive": all_passing_low_positive,
        "affinity": affinity,
        "cpu1_excluded": 1 not in affinity,
        "bootstrap_reps": args.bootstrap_reps,
        "seed": SEED,
        "model_summary": model.to_dict("records"),
        "input_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                D43,
                D104_A,
                D127_OUT / "laysan_last_passage_event_metrics.csv.gz",
                D127_OUT / "joint_model_summary.csv",
                PREREG,
            )
        },
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "One pre-registered summary-statistic sensitivity; median(log CHL) remains the primary bout-level definition.",
    }
    (args.output / "final_summary.json").write_text(
        json.dumps(finite(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"PROGRESS complete 4/4 (100.0%) verdict={payload['verdict']}", flush=True)


if __name__ == "__main__":
    main()
