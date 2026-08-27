#!/usr/bin/env python3
"""Independent event-table and different-seed bootstrap audit of D126."""

from __future__ import annotations

import json
import os
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
D104 = ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz"
D49 = ROOT / "data/audit_inputs/rfbo_classified_event_features.csv.gz"
FORMAL = ROOT / "results/booby_dive_timing/direction_contrasts.csv"
OUT = ROOT / "output/audits/booby_dive_timing"
SCALES = (250, 500, 1000, 2000)
SEED = 126999


def main() -> None:
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 is forbidden; affinity={affinity}")
    OUT.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(D104, usecols=["dataset", "scale_m", "event_id", "cluster", "L_union"], low_memory=False)
    a = a.loc[a.dataset.eq("rfbo") & a.scale_m.isin(SCALES)].copy()
    b = pd.read_csv(D49, usecols=["delta_m", "event_id", "group"], low_memory=False).drop_duplicates()
    joined = a.merge(b, left_on=["scale_m", "event_id"], right_on=["delta_m", "event_id"], how="left", validate="one_to_one")
    if joined.group.isna().any():
        raise RuntimeError("D49 join incomplete")
    formal = pd.read_csv(FORMAL)
    formal = formal.loc[formal.contrast.eq("post_minus_pre") & formal.estimand.eq("L_union")].sort_values("scale_m")
    rng = np.random.default_rng(SEED)
    rows = []
    for scale in SCALES:
        frame = joined.loc[joined.scale_m.eq(scale) & joined.group.isin(["post_only", "pre_only"])]
        bird = frame.groupby(["cluster", "group"]).L_union.mean().unstack()
        paired = bird[["post_only", "pre_only"]].dropna()
        unit = (paired.post_only - paired.pre_only).to_numpy(float)
        observed = float(unit.mean())
        draws = unit[rng.integers(0, len(unit), size=(50000, len(unit)))].mean(axis=1)
        low, high = np.quantile(draws, [0.025, 0.975])
        target = float(formal.loc[formal.scale_m.eq(scale), "observed_unit_equal"].iloc[0])
        rows.append({
            "scale_m": scale, "paired_birds": len(unit), "observed_recomputed": observed,
            "formal_observed": target, "absolute_difference": abs(observed - target),
            "bootstrap_ci_low_seed126999": float(low), "bootstrap_ci_high_seed126999": float(high),
            "ci_crosses_zero": bool(low <= 0 <= high),
        })
    result = pd.DataFrame(rows)
    passed = bool(result.absolute_difference.max() <= 1e-12 and result.ci_crosses_zero.all())
    result.to_csv(OUT / "recomputed_effect_bootstrap.csv", index=False)
    payload = {
        "status": "PASS" if passed else "FAIL", "affinity": affinity, "cpu1_excluded": 1 not in affinity,
        "seed": SEED, "bootstrap_reps": 50000,
        "maximum_observed_absolute_difference": float(result.absolute_difference.max()),
        "all_different_seed_primary_cis_cross_zero": bool(result.ci_crosses_zero.all()),
        "conclusion": "No-temporal-direction verdict is already forced by the bird-bootstrap gate" if passed else "audit discrepancy",
    }
    (OUT / "summary.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False))


if __name__ == "__main__":
    main()
