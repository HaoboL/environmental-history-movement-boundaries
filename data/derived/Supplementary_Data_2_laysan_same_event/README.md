# Supplementary Data 2 laysan same event

This archive contains frozen derived data used in Paper 2. It does not contain rerun experiments.
Empty cells denote quantities not applicable to that row type; they are not silently imputed.

## Files

- `laysan_last_passage_event_metrics.csv.gz` — One row per qualified Laysan event-scale combination. Primary key: `individual_id + event_id + scale_m`.

- `laysan_phase_summary.csv` — One row per Laysan scale-estimand formal phase summary. Primary key: `scale_m + estimand`.

- `laysan_phase_null_audit.csv.gz` — One row per Laysan scale-estimand-phase null estimate. Primary key: `scale_m + estimand + phase`.

- `laysan_reconstruction_audit.csv.gz` — One row per reconstructed Laysan event-scale combination. Primary key: `individual_id + event_id + scale_m`.

- `joint_model_summary.csv` — One row per dataset-scale same-event joint model. Primary key: `dataset + scale_m`.

- `conditional_L_3x3_grid.csv` — One row per dataset-scale-bout-background-tertile-length-tertile cell. Primary key: `dataset + scale_m + abs_tertile + length_tertile`.

- `mean_logchl_joint_model_summary.csv` — One row per dataset-scale joint model after replacing median(log CHL) by mean(log CHL). Primary key: `dataset + scale_m`.

- `mean_logchl_conditional_L_3x3_grid.csv` — One row per dataset-scale-background-tertile-length-tertile cell for the mean(log CHL) sensitivity. Primary key: `dataset + scale_m + abs_tertile + length_tertile`.

- `mean_vs_median_joint_model_comparison.csv` — Side-by-side mean(log CHL) sensitivity and primary median(log CHL) joint-model summaries. Primary key: `dataset + scale_m`.

- `mean_logchl_sensitivity_final_summary.json` — Frozen gate verdict, hashes and model readback for the mean(log CHL) sensitivity. Primary key: `JSON object`.

## Interpretation boundary

CHL is an environmental productivity proxy, not a demonstrated sensory cue. Within-path upper and lower tails are event-relative ranks, not absolute rich and poor habitat. Behavioural dives indicate attempts or underwater activity, not confirmed prey capture.

## Reproducibility

`VARIABLE_DICTIONARY.csv` defines every tabular column. File hashes are listed in `PUBLIC_FILE_MANIFEST.csv` in the parent directory.
