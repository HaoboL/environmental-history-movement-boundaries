# Supplementary Data 1 last passage

This archive contains frozen derived data used in Paper 2. It does not contain rerun experiments.
Empty cells denote quantities not applicable to that row type; they are not silently imputed.

## Files

- `family_a_event_metrics.csv.gz` — One row per qualified dataset-scale-event for E=R+L decomposition. Primary key: `dataset + scale_m + event_id`.

- `family_b_event_metrics.csv.gz` — One row per qualified dataset-scale-event for paired rho-tau localisation. Primary key: `dataset + scale_m + event_id`.

- `formal_summary.csv` — One row per dataset-scale-family-estimand formal summary. Primary key: `dataset + scale_m + family + estimand`.

- `coverage.csv` — Qualification and exclusion counts by dataset and scale. Primary key: `dataset + scale_m + stage`.

## Interpretation boundary

CHL is an environmental productivity proxy, not a demonstrated sensory cue. Within-path upper and lower tails are event-relative ranks, not absolute rich and poor habitat. Behavioural dives indicate attempts or underwater activity, not confirmed prey capture.

## Reproducibility

`VARIABLE_DICTIONARY.csv` defines every tabular column. File hashes are listed in `PUBLIC_FILE_MANIFEST.csv` in the parent directory.
