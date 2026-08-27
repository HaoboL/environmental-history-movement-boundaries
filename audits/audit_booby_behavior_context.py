#!/usr/bin/env python3
"""Independent audit of D106 observed effects and frozen verdict."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
FORMAL = ROOT / "results/booby_behavior_context"
D104 = ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz"
D28 = ROOT / "data/audit_inputs/rfbo_event_two_tail_metrics.csv.gz"
AUDIT_OUT = ROOT / "output/audits/booby_behavior_context"
SCALES = (250, 500, 1000, 2000)
REPS = 20_000
SEED = 1060826


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def holm(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    order = np.argsort(values, kind="stable")
    out = np.ones(len(values), float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(values) - rank) * values[index]))
        out[index] = running
    return out


def boot_ci(unit: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    draws = unit[rng.integers(0, len(unit), size=(REPS, len(unit)))].mean(axis=1)
    return tuple(float(value) for value in np.quantile(draws, [0.025, 0.975]))


def main() -> None:
    AUDIT_OUT.mkdir(parents=True, exist_ok=True)
    affinity = sorted(os.sched_getaffinity(0))
    if 1 in affinity:
        raise RuntimeError(f"CPU1 forbidden in audit: {affinity}")
    d104 = pd.read_csv(D104, low_memory=False)
    d104 = d104.loc[d104.dataset.eq("rfbo") & d104.family.eq("A")].copy()
    d104["event_id"] = d104.event_id.astype(str)
    d28 = pd.read_csv(D28, usecols=["event_id", "delta_m", "behavior_class"], low_memory=False)
    d28["event_id"] = d28.event_id.astype(str)
    d28 = d28.rename(columns={"delta_m": "scale_m"}).drop_duplicates(["scale_m", "event_id"])
    d28 = d28.loc[d28.scale_m.isin(SCALES)].copy()
    merged = d104.merge(d28, on=["scale_m", "event_id"], how="outer", validate="one_to_one", indicator=True)
    if not merged._merge.eq("both").all():
        raise RuntimeError(f"event join failed: {merged._merge.value_counts().to_dict()}")
    identity = float(np.max(np.abs(
        merged[["E_high", "E_low", "E_union"]].to_numpy(float)
        - merged[["R_high", "R_low", "R_union"]].to_numpy(float)
        - merged[["L_high", "L_low", "L_union"]].to_numpy(float)
    )))
    if identity > 1e-12:
        raise RuntimeError("E=R+L audit failed")

    groups = pd.read_csv(FORMAL / "behavior_group_effects.csv")
    modifiers = pd.read_csv(FORMAL / "behavior_modifiers.csv")
    summary = json.loads((FORMAL / "summary.json").read_text(encoding="utf-8"))
    rng = np.random.default_rng(SEED)
    rows = []
    maximum_point_difference = 0.0
    for scale in SCALES:
        frame = merged.loc[merged.scale_m.eq(scale)]
        by_bird = frame.groupby(["cluster", "behavior_class"])["L_union"].mean().unstack()
        dive = by_bird.dive_near.dropna().to_numpy(float)
        paired = by_bird[["dive_near", "wet_only_near"]].dropna()
        difference = (paired.dive_near - paired.wet_only_near).to_numpy(float)
        dive_point = float(dive.mean()); difference_point = float(difference.mean())
        dive_ci = boot_ci(dive, rng); difference_ci = boot_ci(difference, rng)
        output_dive = groups.loc[
            groups.scale_m.eq(scale) & groups.behavior_class.eq("dive_near") & groups.estimand.eq("L_union")
        ].iloc[0]
        output_difference = modifiers.loc[
            modifiers.scale_m.eq(scale) & modifiers.contrast.eq("dive_minus_wet") & modifiers.estimand.eq("L_union")
        ].iloc[0]
        maximum_point_difference = max(
            maximum_point_difference,
            abs(dive_point - float(output_dive.observed_unit_equal)),
            abs(difference_point - float(output_difference.observed_unit_equal)),
        )
        rows.append({
            "scale_m": scale, "dive_birds": len(dive), "paired_birds": len(difference),
            "dive_point": dive_point, "independent_dive_ci_low": dive_ci[0], "independent_dive_ci_high": dive_ci[1],
            "dive_minus_wet_point": difference_point,
            "independent_modifier_ci_low": difference_ci[0], "independent_modifier_ci_high": difference_ci[1],
            "modifier_ci_blocks_primary_pass": difference_ci[0] <= 0,
        })
    audited = pd.DataFrame(rows)
    if maximum_point_difference > 1e-12:
        raise RuntimeError(f"observed point estimate mismatch: {maximum_point_difference}")
    if not audited.modifier_ci_blocks_primary_pass.all():
        raise RuntimeError("independent bootstrap no longer blocks every modifier scale")

    primary_group = groups.loc[groups.behavior_class.eq("dive_near") & groups.estimand.eq("L_union")].sort_values("scale_m")
    primary_modifier = modifiers.loc[modifiers.contrast.eq("dive_minus_wet") & modifiers.estimand.eq("L_union")].sort_values("scale_m")
    group_holm_error = float(np.max(np.abs(holm(primary_group.phase_p_one_sided.to_numpy()) - primary_group.holm_phase_p.to_numpy())))
    modifier_holm_error = float(np.max(np.abs(holm(primary_modifier.phase_p_one_sided.to_numpy()) - primary_modifier.holm_phase_p.to_numpy())))
    if max(group_holm_error, modifier_holm_error) > 1e-12:
        raise RuntimeError("Holm audit failed")
    if primary_modifier.primary_pass.any():
        raise RuntimeError("formal modifier unexpectedly passed")
    if summary["verdict"]["verdict"] != "NO_DIRECT_FEEDING_INTENT_SPECIFICITY":
        raise RuntimeError("summary verdict mismatch")
    if summary["affinity"] != [8] or not summary["cpu1_excluded"]:
        raise RuntimeError("formal CPU affinity audit failed")

    audited.to_csv(AUDIT_OUT / "independent_bootstrap_audit.csv", index=False)
    report = {
        "status": "D106_FORMAL_V1_AUDIT_PASS",
        "independent_seed": SEED,
        "independent_bootstrap_reps": REPS,
        "affinity": affinity,
        "events": int(len(merged)),
        "maximum_E_minus_R_minus_L_error": identity,
        "maximum_observed_point_difference": maximum_point_difference,
        "group_holm_error": group_holm_error,
        "modifier_holm_error": modifier_holm_error,
        "all_four_modifier_CIs_block_primary_pass": True,
        "negative_verdict_does_not_depend_on_phase_p": True,
        "formal_summary_sha256": sha256(FORMAL / "summary.json"),
        "formal_groups_sha256": sha256(FORMAL / "behavior_group_effects.csv"),
        "formal_modifiers_sha256": sha256(FORMAL / "behavior_modifiers.csv"),
        "audit_script_sha256": sha256(Path(__file__)),
    }
    (AUDIT_OUT / "audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("D106_FORMAL_V1_AUDIT_PASS")


if __name__ == "__main__":
    main()
