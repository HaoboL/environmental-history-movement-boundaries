#!/usr/bin/env python3
"""Assemble Paper 2 V3 Source Data and Supplementary Data packages.

The script copies frozen analysis tables and figure-source tables only. It
does not rerun segmentation, environmental extraction, bootstrapping,
randomisation or model fitting.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT
MANUSCRIPT = ROOT / "manuscript"
FIGURES = ROOT / "output/figures"
OUT = ROOT / "output/supplementary_data"

D104 = ROOT / "results/last_record_decomposition"
D127 = ROOT / "results/laysan_same_event"
D45 = ROOT / "results/landing_timing"
D107 = ROOT / "results/shearwater_behavior"
D106 = ROOT / "results/booby_behavior_context"
D126 = ROOT / "results/booby_dive_timing"
D128 = ROOT / "results/boundary_counterfactual"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


PACKAGES: dict[str, list[tuple[Path, str, str, str]]] = {
    "Supplementary_Data_1_last_passage": [
        (D104 / "family_a_event_metrics.csv.gz", "family_a_event_metrics.csv.gz",
         "One row per qualified dataset-scale-event for E=R+L decomposition.",
         "dataset + scale_m + event_id"),
        (D104 / "family_b_event_metrics.csv.gz", "family_b_event_metrics.csv.gz",
         "One row per qualified dataset-scale-event for paired rho-tau localisation.",
         "dataset + scale_m + event_id"),
        (D104 / "formal_summary.csv", "formal_summary.csv",
         "One row per dataset-scale-family-estimand formal summary.",
         "dataset + scale_m + family + estimand"),
        (D104 / "coverage.csv", "coverage.csv",
         "Qualification and exclusion counts by dataset and scale.",
         "dataset + scale_m + stage"),
    ],
    "Supplementary_Data_2_laysan_same_event": [
        (D127 / "laysan_last_passage_event_metrics.csv.gz", "laysan_last_passage_event_metrics.csv.gz",
         "One row per qualified Laysan event-scale combination.",
         "individual_id + event_id + scale_m"),
        (D127 / "laysan_phase_summary.csv", "laysan_phase_summary.csv",
         "One row per Laysan scale-estimand formal phase summary.",
         "scale_m + estimand"),
        (D127 / "laysan_phase_null_audit.csv.gz", "laysan_phase_null_audit.csv.gz",
         "One row per Laysan scale-estimand-phase null estimate.",
         "scale_m + estimand + phase"),
        (D127 / "laysan_reconstruction_audit.csv.gz", "laysan_reconstruction_audit.csv.gz",
         "One row per reconstructed Laysan event-scale combination.",
         "individual_id + event_id + scale_m"),
        (D127 / "joint_model_summary.csv", "joint_model_summary.csv",
         "One row per dataset-scale same-event joint model.",
         "dataset + scale_m"),
        (D127 / "conditional_L_3x3_grid.csv", "conditional_L_3x3_grid.csv",
         "One row per dataset-scale-absolute-CHL-tertile-length-tertile cell.",
         "dataset + scale_m + abs_tertile + length_tertile"),
    ],
    "Supplementary_Data_3_behaviour_context": [
        (D45 / "uesaka_scale_direction_summary.csv", "wandering_albatross_landing_scale_summary.csv",
         "One row per drawdown scale for matched post-minus-pre landing probability.",
         "delta_m"),
        (D45 / "uesaka_lag_bin_summary.csv", "wandering_albatross_landing_timing_summary.csv",
         "One row per drawdown scale and matched temporal bin around a boundary.",
         "delta_m + bin_low_s + bin_high_s"),
        (D107 / "primary_slope_gate.csv", "shearwater_continuous_context_contrast.csv",
         "One row per shearwater scale for continuous foraging-minus-rest slope.",
         "scale_m"),
        (D107 / "behavior_group_effects.csv", "shearwater_behavior_group_effects.csv",
         "One row per shearwater scale-behaviour-group-estimand summary.",
         "scale_m + behavior_group + estimand"),
        (D107 / "event_metrics.csv.gz", "shearwater_event_metrics.csv.gz",
         "One row per qualified shearwater event-scale combination with future behaviour fractions.",
         "individual_id + event_id + scale_m"),
        (D106 / "behavior_group_effects.csv", "red_footed_booby_behavior_group_effects.csv",
         "One row per booby scale-behaviour-class-estimand summary.",
         "scale_m + behavior_class + estimand"),
        (D106 / "behavior_modifiers.csv", "red_footed_booby_behavior_modifiers.csv",
         "One row per booby scale paired dive-versus-wet contrast.",
         "scale_m + contrast + estimand"),
        (D126 / "direction_contrasts.csv", "red_footed_booby_dive_timing_contrasts.csv",
         "One row per booby scale paired post-dive-versus-pre-dive contrast.",
         "scale_m + contrast + estimand"),
    ],
    "Supplementary_Data_4_boundary_counterfactual": [
        (D128 / "rho_segment_metrics.csv", "rho_segment_metrics.csv",
         "One row per observed Crozet segment under the last-record partition.",
         "segment_uid"),
        (D128 / "eligible_record_counterfactual_distribution.csv.gz",
         "eligible_record_counterfactual_distribution.csv.gz",
         "One row per complete eligible-record boundary-replacement replicate.",
         "replicate"),
        (D128 / "segment_boundary_audit.csv", "segment_boundary_audit.csv",
         "One row per segment for observed-partition and confirmation-index diagnostics.",
         "segment_uid"),
        (D128 / "tau_confirmation_reach_diagnostic.csv", "tau_confirmation_reach_diagnostic.csv",
         "Descriptive confirmation-reach statistics; confirmation points are not used as a partition.",
         "reported grouping columns"),
        (D128 / "final_summary.json", "final_summary.json",
         "Observed effects, replacement quantiles and randomisation values.",
         "JSON object"),
    ],
}


SOURCE_DATA: list[tuple[str, Path, str, str]] = [
    ("Fig1", FIGURES / "fig1_source_data.csv", "Main Figure 1, panels a-e",
     "panel + record_type + event_key or family + estimand"),
    ("Fig2", FIGURES / "fig2_source_data.csv", "Main Figure 2, panels a-d",
     "panel + dataset + scale_m + estimand or time bin"),
    ("Fig3", FIGURES / "fig3_source_data.csv", "Main Figure 3, panels a-f",
     "panel + record_type + dataset + scale_m + model term or tertile cell"),
    ("Fig4", FIGURES / "fig4_source_data.csv", "Main Figure 4, panels a-d",
     "panel + record_type + replicate or event position"),
    ("FigS1", FIGURES / "figS1_observation_support_source_data.csv",
     "Supplementary Figure 1, panels a-b", "panel + system + metric + x_value"),
    ("FigS2", FIGURES / "figS2_complete_last_passage_source_data.csv",
     "Supplementary Figure 2, panels a-f", "dataset + scale_m + family + estimand"),
    ("FigS3", FIGURES / "figS3_laysan_conditional_source_data.csv",
     "Supplementary Figure 3, panels a-d", "panel + scale_m + estimand or tertile cell"),
    ("FigS4", FIGURES / "figS4_behaviour_context_source_data.csv",
     "Supplementary Figure 4, panels a-d", "panel + system + scale_m + estimand"),
    ("FigS5", FIGURES / "figS5_boundary_counterfactual_source_data.csv",
     "Supplementary Figure 5, panels a-d", "panel + record_type + replicate or property"),
]


EXACT_DESCRIPTIONS = {
    "dataset": "Short analysis-system identifier.",
    "scale_m": "Radial-drawdown scale in metres.",
    "delta_m": "Radial-drawdown scale in metres.",
    "event_id": "Stable identifier for one movement event within its source trajectory.",
    "event_key": "Stable composite identifier for one movement event.",
    "segment_uid": "Stable identifier for one continuous movement segment.",
    "individual_id": "Stable anonymised biological-individual identifier.",
    "track_id": "Stable source trajectory or deployment identifier.",
    "replicate": "Counterfactual, phase-null or bootstrap replicate identifier.",
    "estimand": "Name of the reported statistical quantity.",
    "family": "Analysis family: decomposition (A) or rho-tau localisation (B).",
    "panel": "Figure panel receiving the row.",
    "record_type": "Semantic row type within a multi-panel source-data table.",
    "events": "Number of qualified event-scale rows contributing to the estimate.",
    "units": "Number of biological sampling units contributing to the estimate.",
    "birds": "Number of birds contributing to the estimate.",
    "individuals": "Number of individuals contributing to the estimate.",
    "observed_unit_equal": "Observed estimate after equal weighting of biological units.",
    "bootstrap_ci_low": "Lower endpoint of the 95% biological-unit-cluster bootstrap interval.",
    "bootstrap_ci_high": "Upper endpoint of the 95% biological-unit-cluster bootstrap interval.",
    "phase_null_mean": "Mean estimate across structure-preserving circular-phase null replicates.",
    "holm_phase_p": "Holm-adjusted one-sided phase-randomisation P value.",
    "logCHL": "Natural logarithm of sampled chlorophyll-a concentration.",
    "Lat": "Latitude in decimal degrees north.",
    "Lon": "Longitude in decimal degrees east.",
    "Time": "Observation timestamp in the source time standard.",
}


def column_description(column: str) -> tuple[str, str]:
    if column in EXACT_DESCRIPTIONS:
        unit = "m" if column in {"scale_m", "delta_m"} else "see description"
        return EXACT_DESCRIPTIONS[column], unit
    lower = column.lower()
    if lower.endswith("_ci_low") or lower.endswith("_ci_high"):
        return "Lower or upper endpoint of the reported 95% confidence interval for the named estimate.", "estimate units"
    if "p" in lower and ("phase" in lower or lower.endswith("_p") or "p_value" in lower):
        return "Randomisation or adjusted probability value for the named comparison.", "probability"
    if "chl" in lower:
        return f"Chlorophyll-a-derived quantity denoted by `{column}`; log fields use natural logarithms and rank fields are within-event fractions.", "field-specific"
    if lower.startswith("l_"):
        return "Last-record component for the tail named in the column suffix.", "probability contrast"
    if lower.startswith("r_"):
        return "Ordinary radial-record component for the tail named in the column suffix.", "probability contrast"
    if lower.startswith("e_"):
        return "Endpoint-excess quantity for the tail named in the column suffix.", "probability contrast"
    if "rho" in lower:
        return f"Quantity evaluated at or involving the retrospectively identified last radial record (`{column}`).", "field-specific"
    if "tau" in lower:
        return f"Quantity evaluated at or involving the later drawdown-confirmation point (`{column}`).", "field-specific"
    if lower.startswith("n_") or lower.endswith("_count") or lower in {"segments", "rows"}:
        return f"Count of {column.replace('_', ' ')}.", "count"
    if "fraction" in lower or "prob" in lower or "percent" in lower:
        return f"Proportion or percentage denoted by {column.replace('_', ' ')}.", "proportion or %"
    if "length" in lower or "dist" in lower:
        return f"Movement-distance or length quantity denoted by {column.replace('_', ' ')}; consult the table README for transformation.", "km, m or log length as named"
    if lower.endswith("_s") or "time" in lower or "lag" in lower:
        return f"Time or lag quantity denoted by {column.replace('_', ' ')}.", "seconds unless stated"
    if "aic" in lower or "support" in lower:
        return f"Model-support quantity denoted by {column.replace('_', ' ')}; positive relative support favours Lomax over exponential.", "AIC units or correlation"
    if "estimate" in lower or "effect" in lower or "contrast" in lower or "slope" in lower or "coefficient" in lower:
        return f"Statistical estimate denoted by {column.replace('_', ' ')}.", "estimate units"
    if lower.startswith("is_") or lower.endswith("_match") or lower.endswith("_valid"):
        return f"Boolean indicator for {column.replace('_', ' ')}.", "0/1 or false/true"
    return f"Frozen analysis field: {column.replace('_', ' ')}.", "field-specific"


def read_table(path: Path) -> pd.DataFrame | None:
    if path.suffix == ".json":
        return None
    return pd.read_csv(path)


def write_variable_dictionary(destination_dir: Path, entries: list[tuple[Path, str, str, str]]) -> None:
    usage: dict[str, dict[str, object]] = defaultdict(lambda: {"files": [], "dtypes": set()})
    for source, filename, _, _ in entries:
        frame = read_table(source)
        if frame is None:
            continue
        for column, dtype in frame.dtypes.items():
            usage[str(column)]["files"].append(filename)
            usage[str(column)]["dtypes"].add(str(dtype))
    rows = []
    for column in sorted(usage):
        description, unit = column_description(column)
        rows.append({
            "column": column,
            "dtype": "; ".join(sorted(usage[column]["dtypes"])),
            "unit": unit,
            "description": description,
            "files": "; ".join(sorted(usage[column]["files"])),
        })
    pd.DataFrame(rows).to_csv(destination_dir / "VARIABLE_DICTIONARY.csv", index=False)


def write_package_readme(destination_dir: Path, package: str,
                         entries: list[tuple[Path, str, str, str]]) -> None:
    lines = [
        f"# {package.replace('_', ' ')}", "",
        "This archive contains frozen derived data used in Paper 2. It does not contain rerun experiments.",
        "Empty cells denote quantities not applicable to that row type; they are not silently imputed.",
        "", "## Files", "",
    ]
    for _, filename, description, key in entries:
        lines.extend([f"- `{filename}` — {description} Primary key: `{key}`.", ""])
    lines.extend([
        "## Interpretation boundary", "",
        "CHL is an environmental productivity proxy, not a demonstrated sensory cue. Within-path upper and lower tails are event-relative ranks, not absolute rich and poor habitat. Behavioural dives indicate attempts or underwater activity, not confirmed prey capture.",
        "", "## Reproducibility", "",
        "`VARIABLE_DICTIONARY.csv` defines every tabular column. File hashes are listed in `PUBLIC_FILE_MANIFEST.csv` in the parent directory.",
    ])
    (destination_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def zip_directory(directory: Path) -> Path:
    archive = OUT / f"{directory.name}.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as handle:
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                handle.write(path, arcname=f"{directory.name}/{path.relative_to(directory)}")
    return archive


def build_source_data() -> list[dict[str, object]]:
    csv_dir = OUT / "source_data_csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    workbook = OUT / "Source_Data.xlsx"
    manifest_rows: list[dict[str, object]] = []
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        for index, (sheet, source, panels, key) in enumerate(SOURCE_DATA, start=1):
            if not source.is_file():
                raise FileNotFoundError(source)
            frame = pd.read_csv(source)
            frame.to_excel(writer, sheet_name=sheet, index=False)
            csv_name = f"{sheet}.csv"
            destination = csv_dir / csv_name
            frame.to_csv(destination, index=False)
            manifest_rows.append({
                "sheet": sheet,
                "csv_file": f"source_data_csv/{csv_name}",
                "figure_panels": panels,
                "row_unit": "One source-data record; panel and record_type identify panel-specific schema.",
                "primary_key": key,
                "rows": len(frame),
                "columns": len(frame.columns),
                "csv_sha256": sha256(destination),
            })
            print(f"[source data {index}/9] {sheet}: {len(frame):,} rows", flush=True)
    pd.DataFrame(manifest_rows).to_csv(OUT / "SOURCE_DATA_MANIFEST.csv", index=False)
    return manifest_rows


def build_description(package_rows: list[dict[str, object]], source_rows: list[dict[str, object]]) -> Path:
    path = MANUSCRIPT / "DESCRIPTION_OF_ADDITIONAL_SUPPLEMENTARY_FILES.md"
    lines = [
        "# Description of additional supplementary files", "",
        "## Environmental history shapes movement boundaries and search scaling", "",
        "The article is accompanied by one Source Data workbook, nine figure-level CSV files and four Supplementary Data archives. All files contain frozen derived results; packaging did not rerun scientific analyses.",
        "", "## Source Data", "",
        "`Source_Data.xlsx` contains one worksheet for each main and supplementary figure. The same values are supplied as CSV files in `source_data_csv`. `SOURCE_DATA_MANIFEST.csv` gives figure mapping, row unit, primary key, dimensions and checksum.",
        "",
    ]
    for row in source_rows:
        lines.append(f"- `{row['sheet']}` — {row['figure_panels']}; {row['rows']:,} rows × {row['columns']} columns; key: `{row['primary_key']}`.")
    lines.extend(["", "## Supplementary Data archives", ""])
    for package in PACKAGES:
        rows = [row for row in package_rows if row["package"] == package]
        lines.extend([
            f"### {package.replace('_', ' ')}", "",
            f"Archive `{package}.zip` contains {len(rows)} frozen analysis tables plus `README.md` and `VARIABLE_DICTIONARY.csv`.", "",
        ])
        for row in rows:
            lines.append(f"- `{row['file']}` — {row['description']} Primary key: `{row['primary_key']}`; {row['rows']} rows when tabular.")
        lines.append("")
    lines.extend([
        "## Integrity and interpretation", "",
        "`PUBLIC_FILE_MANIFEST.csv` reports relative file name, size and SHA-256 checksum without local source paths. `SUPPLEMENTARY_DATA_PACKAGING_AUDIT.json` is an internal provenance record and is not part of the public package. CHL is a productivity proxy rather than a demonstrated cue; path-relative tails are event-relative ranks; behavioural dives do not demonstrate prey capture.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    started = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    print("[0/6 0%] supplementary packaging started", flush=True)
    source_rows = build_source_data()
    print(f"[1/6 17%] source data complete; elapsed={time.monotonic()-started:.1f}s", flush=True)

    package_rows: list[dict[str, object]] = []
    internal_sources: dict[str, dict[str, str]] = {}
    archives: list[Path] = []
    for index, (package, entries) in enumerate(PACKAGES.items(), start=1):
        destination_dir = OUT / package
        destination_dir.mkdir(parents=True, exist_ok=True)
        for source, filename, description, primary_key in entries:
            if not source.is_file():
                raise FileNotFoundError(source)
            destination = destination_dir / filename
            shutil.copy2(source, destination)
            frame = read_table(destination)
            rows = "JSON" if frame is None else int(len(frame))
            package_rows.append({
                "package": package, "file": filename, "description": description,
                "primary_key": primary_key, "rows": rows,
                "bytes": destination.stat().st_size, "sha256": sha256(destination),
            })
            internal_sources[f"{package}/{filename}"] = {
                "source": str(source.relative_to(ROOT)),
                "source_sha256": sha256(source),
                "packaged_sha256": sha256(destination),
            }
        write_variable_dictionary(destination_dir, entries)
        write_package_readme(destination_dir, package, entries)
        archives.append(zip_directory(destination_dir))
        print(f"[{index+1}/6 {17*(index+1)}%] {package} complete; elapsed={time.monotonic()-started:.1f}s", flush=True)

    pd.DataFrame(package_rows).to_csv(OUT / "SUPPLEMENTARY_DATA_CONTENTS.csv", index=False)
    description = build_description(package_rows, source_rows)

    public_paths = [OUT / "Source_Data.xlsx", OUT / "SOURCE_DATA_MANIFEST.csv",
                    OUT / "SUPPLEMENTARY_DATA_CONTENTS.csv", description]
    public_paths.extend(sorted((OUT / "source_data_csv").glob("*.csv")))
    public_paths.extend(archives)
    public_rows = [{"file": str(path.relative_to(MANUSCRIPT)), "bytes": path.stat().st_size,
                    "sha256": sha256(path)} for path in public_paths]
    public_manifest_path = OUT / "PUBLIC_FILE_MANIFEST.csv"
    with public_manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(public_rows)

    audit = {
        "status": "complete", "source_data_sheets": len(source_rows),
        "supplementary_archives": len(archives),
        "packaged_analysis_files": len(package_rows),
        "scientific_recomputation": False, "cpu1_used": False,
        "elapsed_seconds": time.monotonic() - started,
        "inputs": internal_sources,
        "public_manifest_sha256": sha256(public_manifest_path),
    }
    (OUT / "SUPPLEMENTARY_DATA_PACKAGING_AUDIT.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[6/6 100%] packaging complete; elapsed={time.monotonic()-started:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
