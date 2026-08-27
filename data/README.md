# Data scope and provenance

This repository contains derived, article-specific tables—not the original animal observations or complete satellite products.

`source_data/` contains the CSV equivalents of the data plotted in Figs. 1–4 and Supplementary Figs. 1–5. `derived/` contains the submitted `Source_Data.xlsx`, variable dictionaries, exclusion/coverage information and Supplementary Data tables. `audit_inputs/` contains two compact derived event tables needed to rerun independent behaviour audits. `results/` at repository root contains the frozen formal output objects read by the release checks and figure builders.

The original public datasets and their persistent identifiers are listed in the top-level README. A full raw-to-result rerun requires downloading those records from their source repositories and placing the resulting inputs in the layout referenced by the portable analysis scripts. No API keys, account credentials, private animal data or proprietary vessel data are included.

Derived tables retain source event or deployment identifiers where those identifiers are necessary to audit biological-unit clustering and one-to-one joins. They do not add confidential locations beyond the public source records.
