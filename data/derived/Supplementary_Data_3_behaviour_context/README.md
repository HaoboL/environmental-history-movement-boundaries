# Supplementary Data 3 behaviour context

This archive contains frozen derived data used in Paper 2. It does not contain rerun experiments.
Empty cells denote quantities not applicable to that row type; they are not silently imputed.

## Files

- `wandering_albatross_landing_scale_summary.csv` — One row per drawdown scale for matched post-minus-pre landing probability. Primary key: `delta_m`.

- `wandering_albatross_landing_timing_summary.csv` — One row per drawdown scale and matched temporal bin around a boundary. Primary key: `delta_m + bin_low_s + bin_high_s`.

- `shearwater_continuous_context_contrast.csv` — One row per shearwater scale for continuous foraging-minus-rest slope. Primary key: `scale_m`.

- `shearwater_behavior_group_effects.csv` — One row per shearwater scale-behaviour-group-estimand summary. Primary key: `scale_m + behavior_group + estimand`.

- `shearwater_event_metrics.csv.gz` — One row per qualified shearwater event-scale combination with future behaviour fractions. Primary key: `individual_id + event_id + scale_m`.

- `red_footed_booby_behavior_group_effects.csv` — One row per booby scale-behaviour-class-estimand summary. Primary key: `scale_m + behavior_class + estimand`.

- `red_footed_booby_behavior_modifiers.csv` — One row per booby scale paired dive-versus-wet contrast. Primary key: `scale_m + contrast + estimand`.

- `red_footed_booby_dive_timing_contrasts.csv` — One row per booby scale paired post-dive-versus-pre-dive contrast. Primary key: `scale_m + contrast + estimand`.

## Interpretation boundary

CHL is an environmental productivity proxy, not a demonstrated sensory cue. Within-path upper and lower tails are event-relative ranks, not absolute rich and poor habitat. Behavioural dives indicate attempts or underwater activity, not confirmed prey capture.

## Reproducibility

`VARIABLE_DICTIONARY.csv` defines every tabular column. File hashes are listed in `PUBLIC_FILE_MANIFEST.csv` in the parent directory.
