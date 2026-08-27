# D107 short-tailed shearwater last-passage CHL × behavior result

- Verdict: `NO_DIRECT_FORAGING_STATE_LINK`
- Eligible events: 254; individuals: 9
- Primary slope pass scales: []
- Primary forage-group pass scales: [500, 1000, 2000]
- Adjacent joint-pass pairs: []

## Primary continuous interaction

|   scale_m | predictor          | estimand   |   events |   individuals |   individual_equal_slope |   bootstrap_ci_low |   bootstrap_ci_high |   phase_null_mean |   phase_p_raw | holm_family                                 |   phase_p_holm | coverage_pass   | pass   |
|----------:|:-------------------|:-----------|---------:|--------------:|-------------------------:|-------------------:|--------------------:|------------------:|--------------:|:--------------------------------------------|---------------:|:----------------|:-------|
|       500 | future_specificity | L_union    |       74 |             7 |                0.261099  |         -0.011293  |            0.558095 |        0.0010497  |         0.042 | primary_future_specificity_L_union_4_scales |          0.136 | True            | False  |
|      1000 | future_specificity | L_union    |       70 |             7 |                0.289271  |         -0.0221191 |            0.616708 |       -0.00343979 |         0.034 | primary_future_specificity_L_union_4_scales |          0.136 | True            | False  |
|      2000 | future_specificity | L_union    |       62 |             7 |                0.0978082 |         -0.229866  |            0.366536 |        0.00382543 |         0.313 | primary_future_specificity_L_union_4_scales |          0.626 | True            | False  |
|      5000 | future_specificity | L_union    |       35 |             6 |                0.0880413 |         -0.492966  |            0.959996 |       -0.0501198  |         0.424 | primary_future_specificity_L_union_4_scales |          0.626 | True            | False  |

## Primary forage-dominant within-group footprint

|   scale_m | behavior_group   | estimand   |   events |   individuals |   individual_equal_mean |   bootstrap_ci_low |   bootstrap_ci_high |   phase_null_mean |   phase_p_raw | holm_family                              |   phase_p_holm | coverage_pass   | pass   |
|----------:|:-----------------|:-----------|---------:|--------------:|------------------------:|-------------------:|--------------------:|------------------:|--------------:|:-----------------------------------------|---------------:|:----------------|:-------|
|       500 | forage_dominant  | L_union    |       35 |             9 |                0.356368 |          0.212184  |            0.500424 |         0.0803169 |         0.003 | primary_forage_dominant_L_union_4_scales |          0.012 | True            | True   |
|      1000 | forage_dominant  | L_union    |       34 |             8 |                0.334976 |          0.187709  |            0.484875 |         0.0506222 |         0.005 | primary_forage_dominant_L_union_4_scales |          0.015 | True            | True   |
|      2000 | forage_dominant  | L_union    |       32 |             8 |                0.293344 |          0.115788  |            0.471378 |         0.0449848 |         0.009 | primary_forage_dominant_L_union_4_scales |          0.018 | True            | True   |
|      5000 | forage_dominant  | L_union    |       18 |             8 |                0.19894  |         -0.0032828 |            0.397803 |         0.0462059 |         0.118 | primary_forage_dominant_L_union_4_scales |          0.118 | True            | False  |

## Interpretation boundary

A positive gate supports a same-event link between the CHL last-passage footprint and foraging-versus-rest state. It does not establish direct CHL perception, a sensory channel, successful prey capture, or causal triggering.
