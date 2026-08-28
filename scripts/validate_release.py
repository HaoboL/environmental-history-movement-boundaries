#!/usr/bin/env python3
"""Validate public-release completeness and article-level algebraic invariants."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata/FILE_MANIFEST.csv"
EXCLUDED_PARTS = {
    ".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".venv",
    "output", "external_inputs",
}
TEXT_SUFFIXES = {".py", ".md", ".txt", ".toml", ".yml", ".yaml", ".cff"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def iter_release_files():
    """Yield release files without descending into excluded directory trees."""
    for directory, subdirectories, filenames in os.walk(ROOT):
        subdirectories[:] = [name for name in subdirectories if name not in EXCLUDED_PARTS]
        base = Path(directory)
        for filename in filenames:
            path = base / filename
            if path != MANIFEST:
                yield path


def release_files() -> set[str]:
    return {path.relative_to(ROOT).as_posix() for path in iter_release_files()}


def validate_manifest() -> None:
    table = pd.read_csv(MANIFEST)
    if table.path.duplicated().any():
        raise AssertionError("duplicate path in FILE_MANIFEST.csv")
    expected = set(table.path.astype(str))
    actual = release_files()
    if expected != actual:
        raise AssertionError(
            f"manifest membership mismatch; missing={sorted(expected-actual)}, extra={sorted(actual-expected)}"
        )
    for row in table.itertuples(index=False):
        path = ROOT / row.path
        if path.stat().st_size != int(row.bytes):
            raise AssertionError(f"byte-count mismatch: {row.path}")
        if digest(path) != row.sha256:
            raise AssertionError(f"SHA-256 mismatch: {row.path}")


def validate_last_record_identities() -> None:
    a = pd.read_csv(ROOT / "results/last_record_decomposition/family_a_event_metrics.csv.gz")
    b = pd.read_csv(ROOT / "results/last_record_decomposition/family_b_event_metrics.csv.gz")
    a_error = np.abs(
        a[["E_high", "E_low", "E_union"]].to_numpy(float)
        - a[["R_high", "R_low", "R_union"]].to_numpy(float)
        - a[["L_high", "L_low", "L_union"]].to_numpy(float)
    ).max()
    b_error = np.abs(
        b[["C_high", "C_low", "C_union"]].to_numpy(float)
        - b[["rho_high", "rho_low", "rho_union"]].to_numpy(float)
        + b[["tau_high", "tau_low", "tau_union"]].to_numpy(float)
    ).max()
    if a_error > 1e-12 or b_error > 1e-12:
        raise AssertionError(f"algebraic identity failed: E=R+L {a_error}; C=rho-tau {b_error}")
    if a.duplicated(["dataset", "scale_m", "event_id"]).any():
        raise AssertionError("duplicate Family-A event key")
    if b.duplicated(["dataset", "scale_m", "event_id"]).any():
        raise AssertionError("duplicate Family-B event key")


def validate_source_data() -> None:
    source_dir = ROOT / "data/source_data"
    expected = {f"Fig{i}.csv" for i in range(1, 5)} | {f"FigS{i}.csv" for i in range(1, 6)}
    actual = {path.name for path in source_dir.glob("*.csv")}
    aliases = {
        "fig1_source_data.csv": "Fig1.csv", "fig2_source_data.csv": "Fig2.csv",
        "fig3_source_data.csv": "Fig3.csv", "fig4_source_data.csv": "Fig4.csv",
        "figS1_observation_support_source_data.csv": "FigS1.csv",
        "figS2_complete_last_passage_source_data.csv": "FigS2.csv",
        "figS3_laysan_conditional_source_data.csv": "FigS3.csv",
        "figS4_behaviour_context_source_data.csv": "FigS4.csv",
        "figS5_boundary_counterfactual_source_data.csv": "FigS5.csv",
    }
    mapped = {aliases.get(name, name) for name in actual}
    if not expected.issubset(mapped):
        raise AssertionError(f"missing figure source data: {sorted(expected-mapped)}")
    workbook = ROOT / "data/derived/Source_Data.xlsx"
    if not workbook.is_file() or workbook.stat().st_size == 0:
        raise AssertionError("Source_Data.xlsx is absent or empty")
    with zipfile.ZipFile(workbook) as archive:
        workbook_xml = ET.fromstring(archive.read("xl/workbook.xml"))
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    sheet_names = {
        sheet.attrib["name"]
        for sheet in workbook_xml.findall("x:sheets/x:sheet", namespace)
    }
    required_sheets = {Path(name).stem for name in expected}
    missing = required_sheets - sheet_names
    if missing:
        raise AssertionError(f"Source_Data.xlsx missing sheets: {sorted(missing)}")


def validate_mean_logchl_sensitivity() -> None:
    result_dir = ROOT / "results/mean_logchl_sensitivity"
    summary = pd.read_csv(result_dir / "mean_logchl_joint_model_summary.csv")
    if len(summary) != 4:
        raise AssertionError(f"mean(log CHL) sensitivity expected 4 rows, found {len(summary)}")
    required = {
        ("goto", 100),
        ("usgs_laysan_albatross", 500),
        ("usgs_laysan_albatross", 1000),
        ("usgs_laysan_albatross", 2000),
    }
    observed = set(zip(summary.dataset.astype(str), summary.scale_m.astype(int)))
    if observed != required:
        raise AssertionError(f"mean(log CHL) sensitivity systems/scales differ: {sorted(observed)}")
    gates = (
        summary.beta_absolute_ci_high.lt(0)
        & summary.standardized_L_union_ci_low.gt(0)
        & summary.beta_L_low.gt(0)
        & summary.positive_grids.eq(9)
        & summary.estimable_grids.eq(9)
    )
    if not gates.all():
        raise AssertionError("mean(log CHL) sensitivity frozen gate is not satisfied in every row")
    final = json.loads((result_dir / "final_summary.json").read_text(encoding="utf-8"))
    if final.get("verdict") != "MEAN_LOGCHL_SENSITIVITY_PASSED" or not final.get("gate_pass"):
        raise AssertionError("mean(log CHL) sensitivity final verdict is not PASS")
    if int(final.get("bootstrap_reps", 0)) != 20000 or not final.get("cpu1_excluded"):
        raise AssertionError("mean(log CHL) sensitivity resource/readback mismatch")
    archive = ROOT / "data/derived/Supplementary_Data_2_laysan_same_event.zip"
    required_members = {
        "mean_logchl_joint_model_summary.csv",
        "mean_logchl_conditional_L_3x3_grid.csv",
        "mean_vs_median_joint_model_comparison.csv",
        "mean_logchl_sensitivity_final_summary.json",
    }
    with zipfile.ZipFile(archive) as bundle:
        names = {Path(name).name for name in bundle.namelist()}
    if not required_members.issubset(names):
        raise AssertionError(f"Supplementary Data 2 missing mean-sensitivity files: {sorted(required_members-names)}")


def validate_text_safety() -> None:
    forbidden = [
        re.compile("/home/" + "liu/"),
        re.compile(r"(?i)(api[_-]?key|password|access[_-]?token)\s*[:=]\s*['\"]?[A-Za-z0-9_-]{12,}"),
        re.compile(r"github_pat_[A-Za-z0-9_]+"),
        re.compile(r"ghp_[A-Za-z0-9]+"),
    ]
    for path in iter_release_files():
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if "frozen_originals" in path.parts:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for pattern in forbidden:
            if pattern.search(text):
                raise AssertionError(f"potential private path or secret pattern in {relative}")


def main() -> int:
    checks = [
        ("manifest", validate_manifest),
        ("last-record identities", validate_last_record_identities),
        ("source-data coverage", validate_source_data),
        ("mean-logCHL sensitivity", validate_mean_logchl_sensitivity),
        ("text safety", validate_text_safety),
    ]
    for index, (name, check) in enumerate(checks, start=1):
        check()
        print(f"[{index}/{len(checks)}] PASS: {name}")
    print("RELEASE_VALIDATION_PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"RELEASE_VALIDATION_FAIL: {error}", file=sys.stderr)
        raise
