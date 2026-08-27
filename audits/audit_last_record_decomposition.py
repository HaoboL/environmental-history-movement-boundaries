#!/usr/bin/env python3
"""Independent integrity audit for the frozen D104 formal-v2 output.

This script does not rerun RD, bootstrap, or phase randomization.  It checks the
completed artifacts against the preregistered algebra, test families, frozen
event universe, hashes, and decision rule.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER2_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
OUT = ROOT / "results/last_record_decomposition"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def adjacent_scales(scales: set[int]) -> bool:
    ordered = [250, 500, 1000, 2000]
    return any(a in scales and b in scales for a, b in zip(ordered, ordered[1:]))


def main() -> None:
    final = json.loads((OUT / "final_summary.json").read_text())
    summary = pd.read_csv(OUT / "formal_summary.csv")
    coverage = pd.read_csv(OUT / "coverage.csv")
    audit = pd.read_csv(OUT / "event_structure_audit.csv.gz")
    fam_a = pd.read_csv(OUT / "family_a_event_metrics.csv.gz")
    fam_b = pd.read_csv(OUT / "family_b_event_metrics.csv.gz")
    phase = pd.read_csv(OUT / "phase_null_audit.csv.gz")

    require(final["status"] == "FORMAL_COMPLETE", "formal run is not complete")
    require(final["affinity"] == [8], f"unexpected CPU affinity: {final['affinity']}")
    require(final["cpu1_excluded"] is True, "CPU1 exclusion flag is false")
    require(final["phase_reps"] == 999, "phase replicate count is not 999")
    require(final["bootstrap_reps"] == 5000, "bootstrap replicate count is not 5000")
    require(len(summary) == 90, f"expected 90 summary rows, got {len(summary)}")
    require(not summary.duplicated(["dataset", "scale_m", "family", "estimand"]).any(), "duplicate summary rows")
    require(not fam_a.duplicated(["dataset", "scale_m", "event_id"]).any(), "duplicate Family-A events")
    require(not fam_b.duplicated(["dataset", "scale_m", "event_id"]).any(), "duplicate Family-B events")

    expected_rows = {("goto", "A"): 9, ("goto", "B"): 9, ("rfbo", "A"): 36, ("rfbo", "B"): 36}
    got_rows = summary.groupby(["dataset", "family"]).size().to_dict()
    require(got_rows == expected_rows, f"unexpected summary family sizes: {got_rows}")
    expected_tests = {("goto", "A"): 6, ("goto", "B"): 9, ("rfbo", "A"): 24, ("rfbo", "B"): 36}
    got_tests = summary[summary.holm_phase_p.notna()].groupby(["dataset", "family"]).size().to_dict()
    require(got_tests == expected_tests, f"unexpected Holm family sizes: {got_tests}")

    a_err = float(np.abs(fam_a[["E_high", "E_low", "E_union"]].to_numpy()
                         - fam_a[["R_high", "R_low", "R_union"]].to_numpy()
                         - fam_a[["L_high", "L_low", "L_union"]].to_numpy()).max())
    b_err = float(np.abs(fam_b[["C_high", "C_low", "C_union"]].to_numpy()
                         - fam_b[["rho_high", "rho_low", "rho_union"]].to_numpy()
                         + fam_b[["tau_high", "tau_low", "tau_union"]].to_numpy()).max())
    require(a_err <= 1e-12, f"Family-A identity failed: {a_err}")
    require(b_err <= 1e-12, f"Family-B identity failed: {b_err}")
    require(abs(a_err - final["family_a_identity_max_abs_error"]) <= 1e-15, "Family-A reported error mismatch")
    require(abs(b_err - final["family_b_identity_max_abs_error"]) <= 1e-15, "Family-B reported error mismatch")

    for dataset in ("goto", "rfbo"):
        anchor = final["anchor_reproduction"][dataset]
        require(anchor["rows"] == anchor["old_rows_matched"], f"{dataset} anchor row mismatch")
        require(anchor["exact_within_1e_12_fraction"] == 1.0, f"{dataset} anchor not exact")
        require(anchor["maximum_absolute_E_difference"] <= 1e-12, f"{dataset} anchor numerical mismatch")

    require((coverage.geometry_pass == coverage.frozen_events).all(), "not all frozen events passed geometry")
    require((coverage.rho_tau_distinct_cell_fraction.between(0, 1)).all(), "invalid distinct-cell fractions")
    rfbo_audit = audit[(audit.dataset == "rfbo") & audit.a_decomposition_eligible]
    require(rfbo_audit.a_in_frozen_d28_universe.fillna(False).all(), "RFBO Family-A includes events outside D28 universe")
    require((phase.valid_unique_offsets > 0).all(), "a phase-null segment has no valid nonzero offset")
    require((phase.valid_unique_offsets <= phase.attempted_unique_offsets).all(), "valid phase offsets exceed attempts")

    pscaled = summary.phase_p_one_sided.to_numpy() * 1000
    require(np.allclose(pscaled, np.round(pscaled), atol=1e-10), "phase p-values are not on the 1/1000 grid")
    require((summary.phase_p_one_sided >= 0.001).all(), "phase p-value below attainable minimum")
    tested = summary.holm_phase_p.notna()
    expected_support = tested & (summary.bootstrap_ci_low > 0) & (summary.holm_phase_p <= 0.05)
    require((summary.positive_support.astype(bool) == expected_support).all(), "positive-support flags violate frozen rule")

    missing_external_inputs = []
    for relative, expected_hash in final["input_sha256"].items():
        path = ROOT / relative
        if not path.is_file():
            missing_external_inputs.append(relative)
            continue
        require(sha256(path) == expected_hash, f"input hash mismatch: {relative}")

    json_summary = pd.DataFrame(final["summary"])[summary.columns]
    pd.testing.assert_frame_equal(summary, json_summary, check_dtype=False, atol=1e-15, rtol=1e-12)
    json_coverage = pd.DataFrame(final["coverage"])[coverage.columns]
    pd.testing.assert_frame_equal(coverage, json_coverage, check_dtype=False, atol=1e-15, rtol=1e-12)

    # Recompute only the frozen cross-system decisions used in final_summary.
    a_recomputed = {}
    for tail in ("high", "low", "union"):
        components = {}
        for component in ("R", "L"):
            name = f"{component}_{tail}"
            goto_pass = bool(summary[(summary.dataset == "goto") & (summary.family == "A") &
                                     (summary.estimand == name)].positive_support.iloc[0])
            rfbo_scales = set(summary[(summary.dataset == "rfbo") & (summary.family == "A") &
                                      (summary.estimand == name) & summary.positive_support].scale_m.astype(int))
            components[component] = goto_pass and adjacent_scales(rfbo_scales)
        a_recomputed[tail] = components
        require(components == final["verdict"]["family_a"][tail]["component_cross_system_pass"],
                f"Family-A verdict mismatch for {tail}")

    for tail in ("high", "low", "union"):
        name = f"C_{tail}"
        goto_pass = bool(summary[(summary.dataset == "goto") & (summary.family == "B") &
                                 (summary.estimand == name)].positive_support.iloc[0])
        rfbo_scales = sorted(summary[(summary.dataset == "rfbo") & (summary.family == "B") &
                                     (summary.estimand == name) & summary.positive_support].scale_m.astype(int).tolist())
        require(goto_pass == final["verdict"]["family_b"][tail]["goto_pass"], f"Family-B Goto verdict mismatch: {tail}")
        require(rfbo_scales == final["verdict"]["family_b"][tail]["rfbo_pass_scales"],
                f"Family-B RFBO verdict mismatch: {tail}")

    print("LAST_RECORD_DECOMPOSITION_AUDIT_PASS")
    print(f"summary_rows={len(summary)} family_a_identity_max={a_err:.3e} family_b_identity_max={b_err:.3e}")
    print("affinity=8 cpu1_excluded=true phase_reps=999 bootstrap_reps=5000")
    input_status = "PASS" if not missing_external_inputs else f"NOT_RUN({len(missing_external_inputs)} external inputs not redistributed)"
    print(f"input_hashes={input_status} anchors=PASS frozen_verdicts=PASS")


if __name__ == "__main__":
    main()
