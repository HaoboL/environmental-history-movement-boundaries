# Description of additional supplementary files

## Environmental history shapes movement boundaries and search scaling

The article is accompanied by one Source Data workbook, nine figure-level CSV files and four Supplementary Data archives. All files contain frozen derived results; packaging did not rerun scientific analyses.

## Source Data

`Source_Data.xlsx` contains one worksheet for each main and supplementary figure. The same values are supplied as CSV files in `source_data_csv`. `SOURCE_DATA_MANIFEST.csv` gives figure mapping, row unit, primary key, dimensions and checksum.

- `Fig1` — Main Figure 1, panels a-e; 40 rows × 42 columns; key: `panel + record_type + event_key or family + estimand`.
- `Fig2` — Main Figure 2, panels a-d; 40 rows × 33 columns; key: `panel + dataset + scale_m + estimand or time bin`.
- `Fig3` — Main Figure 3, panels a-f; 52 rows × 16 columns; key: `panel + record_type + dataset + scale_m + model term or tertile cell`.
- `Fig4` — Main Figure 4, panels a-d; 1,022 rows × 45 columns; key: `panel + record_type + replicate or event position`.
- `FigS1` — Supplementary Figure 1, panels a-b; 1,703 rows × 7 columns; key: `panel + system + metric + x_value`.
- `FigS2` — Supplementary Figure 2, panels a-f; 90 rows × 17 columns; key: `dataset + scale_m + family + estimand`.
- `FigS3` — Supplementary Figure 3, panels a-d; 33 rows × 41 columns; key: `panel + scale_m + estimand or tertile cell`.
- `FigS4` — Supplementary Figure 4, panels a-d; 35 rows × 12 columns; key: `panel + system + scale_m + estimand`.
- `FigS5` — Supplementary Figure 5, panels a-d; 1,003 rows × 10 columns; key: `panel + record_type + replicate or property`.

## Supplementary Data archives

### Supplementary Data 1 last passage

Archive `Supplementary_Data_1_last_passage.zip` contains 4 frozen analysis tables plus `README.md` and `VARIABLE_DICTIONARY.csv`.

- `family_a_event_metrics.csv.gz` — One row per qualified dataset-scale-event for E=R+L decomposition. Primary key: `dataset + scale_m + event_id`; 16865 rows when tabular.
- `family_b_event_metrics.csv.gz` — One row per qualified dataset-scale-event for paired rho-tau localisation. Primary key: `dataset + scale_m + event_id`; 18689 rows when tabular.
- `formal_summary.csv` — One row per dataset-scale-family-estimand formal summary. Primary key: `dataset + scale_m + family + estimand`; 90 rows when tabular.
- `coverage.csv` — Qualification and exclusion counts by dataset and scale. Primary key: `dataset + scale_m + stage`; 5 rows when tabular.

### Supplementary Data 2 laysan same event

Archive `Supplementary_Data_2_laysan_same_event.zip` contains 10 frozen analysis tables plus `README.md` and `VARIABLE_DICTIONARY.csv`.

- `laysan_last_passage_event_metrics.csv.gz` — One row per qualified Laysan event-scale combination. Primary key: `individual_id + event_id + scale_m`; 10572 rows when tabular.
- `laysan_phase_summary.csv` — One row per Laysan scale-estimand formal phase summary. Primary key: `scale_m + estimand`; 27 rows when tabular.
- `laysan_phase_null_audit.csv.gz` — One row per Laysan scale-estimand-phase null estimate. Primary key: `scale_m + estimand + phase`; 3930 rows when tabular.
- `laysan_reconstruction_audit.csv.gz` — One row per reconstructed Laysan event-scale combination. Primary key: `individual_id + event_id + scale_m`; 10572 rows when tabular.
- `joint_model_summary.csv` — One row per dataset-scale same-event joint model. Primary key: `dataset + scale_m`; 4 rows when tabular.
- `conditional_L_3x3_grid.csv` — One row per dataset-scale-bout-background-tertile-length-tertile cell. Primary key: `dataset + scale_m + abs_tertile + length_tertile`; 36 rows when tabular.
- `mean_logchl_joint_model_summary.csv` — One row per dataset-scale joint model after replacing median(log CHL) by mean(log CHL). Primary key: `dataset + scale_m`; 4 rows when tabular.
- `mean_logchl_conditional_L_3x3_grid.csv` — One row per dataset-scale-background-tertile-length-tertile cell for the mean(log CHL) sensitivity. Primary key: `dataset + scale_m + abs_tertile + length_tertile`; 36 rows when tabular.
- `mean_vs_median_joint_model_comparison.csv` — Side-by-side mean(log CHL) sensitivity and primary median(log CHL) joint-model summaries. Primary key: `dataset + scale_m`; 4 rows when tabular.
- `mean_logchl_sensitivity_final_summary.json` — Frozen gate verdict, hashes and model readback for the mean(log CHL) sensitivity. Primary key: `JSON object`; JSON rows when tabular.

### Supplementary Data 3 behaviour context

Archive `Supplementary_Data_3_behaviour_context.zip` contains 8 frozen analysis tables plus `README.md` and `VARIABLE_DICTIONARY.csv`.

- `wandering_albatross_landing_scale_summary.csv` — One row per drawdown scale for matched post-minus-pre landing probability. Primary key: `delta_m`; 15 rows when tabular.
- `wandering_albatross_landing_timing_summary.csv` — One row per drawdown scale and matched temporal bin around a boundary. Primary key: `delta_m + bin_low_s + bin_high_s`; 60 rows when tabular.
- `shearwater_continuous_context_contrast.csv` — One row per shearwater scale for continuous foraging-minus-rest slope. Primary key: `scale_m`; 4 rows when tabular.
- `shearwater_behavior_group_effects.csv` — One row per shearwater scale-behaviour-group-estimand summary. Primary key: `scale_m + behavior_group + estimand`; 72 rows when tabular.
- `shearwater_event_metrics.csv.gz` — One row per qualified shearwater event-scale combination with future behaviour fractions. Primary key: `individual_id + event_id + scale_m`; 254 rows when tabular.
- `red_footed_booby_behavior_group_effects.csv` — One row per booby scale-behaviour-class-estimand summary. Primary key: `scale_m + behavior_class + estimand`; 108 rows when tabular.
- `red_footed_booby_behavior_modifiers.csv` — One row per booby scale paired dive-versus-wet contrast. Primary key: `scale_m + contrast + estimand`; 72 rows when tabular.
- `red_footed_booby_dive_timing_contrasts.csv` — One row per booby scale paired post-dive-versus-pre-dive contrast. Primary key: `scale_m + contrast + estimand`; 108 rows when tabular.

### Supplementary Data 4 boundary counterfactual

Archive `Supplementary_Data_4_boundary_counterfactual.zip` contains 5 frozen analysis tables plus `README.md` and `VARIABLE_DICTIONARY.csv`.

- `rho_segment_metrics.csv` — One row per observed Crozet segment under the last-record partition. Primary key: `segment_uid`; 129 rows when tabular.
- `eligible_record_counterfactual_distribution.csv.gz` — One row per complete eligible-record boundary-replacement replicate. Primary key: `replicate`; 999 rows when tabular.
- `segment_boundary_audit.csv` — One row per segment for observed-partition and confirmation-index diagnostics. Primary key: `segment_uid`; 129 rows when tabular.
- `tau_confirmation_reach_diagnostic.csv` — Descriptive confirmation-reach statistics; confirmation points are not used as a partition. Primary key: `reported grouping columns`; 129 rows when tabular.
- `final_summary.json` — Observed effects, replacement quantiles and randomisation values. Primary key: `JSON object`; JSON rows when tabular.

## Integrity and interpretation

`PUBLIC_FILE_MANIFEST.csv` reports relative file name, size and SHA-256 checksum without local source paths. `SUPPLEMENTARY_DATA_PACKAGING_AUDIT.json` is an internal provenance record and is not part of the public package. CHL is a productivity proxy rather than a demonstrated cue; path-relative tails are event-relative ranks; behavioural dives do not demonstrate prey capture.
