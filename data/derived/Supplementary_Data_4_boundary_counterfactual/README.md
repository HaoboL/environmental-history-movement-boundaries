# Supplementary Data 4 boundary counterfactual

This archive contains frozen derived data used in Paper 2. It does not contain rerun experiments.
Empty cells denote quantities not applicable to that row type; they are not silently imputed.

## Files

- `rho_segment_metrics.csv` — One row per observed Crozet segment under the last-record partition. Primary key: `segment_uid`.

- `eligible_record_counterfactual_distribution.csv.gz` — One row per complete eligible-record boundary-replacement replicate. Primary key: `replicate`.

- `segment_boundary_audit.csv` — One row per segment for observed-partition and confirmation-index diagnostics. Primary key: `segment_uid`.

- `tau_confirmation_reach_diagnostic.csv` — Descriptive confirmation-reach statistics; confirmation points are not used as a partition. Primary key: `reported grouping columns`.

- `final_summary.json` — Observed effects, replacement quantiles and randomisation values. Primary key: `JSON object`.

## Interpretation boundary

CHL is an environmental productivity proxy, not a demonstrated sensory cue. Within-path upper and lower tails are event-relative ranks, not absolute rich and poor habitat. Behavioural dives indicate attempts or underwater activity, not confirmed prey capture.

## Reproducibility

`VARIABLE_DICTIONARY.csv` defines every tabular column. File hashes are listed in `PUBLIC_FILE_MANIFEST.csv` in the parent directory.
