# Data scope and provenance

This repository contains derived, article-specific tables—not the original animal observations or complete satellite products.

`source_data/` contains the CSV equivalents of the data plotted in Figs. 1–4 and Supplementary Figs. 1–5. `derived/` contains the submitted `Source_Data.xlsx`, variable dictionaries, exclusion/coverage information and Supplementary Data tables. `audit_inputs/` contains two compact derived event tables needed to rerun independent behaviour audits. `results/` at repository root contains the frozen formal output objects read by the release checks and figure builders.

The original public datasets, persistent identifiers, transfer sizes and executable download commands are listed in the top-level README. `scripts/download_public_data.py` retrieves the public animal records, while `scripts/download_chl.py` executes the frozen track-specific Copernicus request manifests. Both write beneath the ignored `external_inputs/` tree so raw observations and credentials cannot be committed accidentally. No API keys, account credentials, private animal data or proprietary vessel data are included.

The Git release and the larger submission archive serve different roles. Git contains the executable code, frozen formal results and every plotted value. The submission-matched permanent archive will additionally freeze the larger event-level analysis inputs needed to rerun formal models without repeating source-specific telemetry cleaning. This keeps the event population used for inference immutable while retaining a documented route back to every original public record.

Derived tables retain source event or deployment identifiers where those identifiers are necessary to audit biological-unit clustering and one-to-one joins. They do not add confidential locations beyond the public source records.
