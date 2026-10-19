# Two-way ZIP-and-season bootstrap for Figure 1

The main point estimates, corrected 0–5 suppression definition, suppressed-group mean of 3.8, RMW swaths, phase calendar, population adjustment and prediction target are retained. Only the resampling design changes.

Within each threshold, 999 bootstrap draws independently resampled all ZIPs and the 18 calendar-year hurricane seasons with replacement. A cell's frequency weight was multiplied by its ZIP multiplicity and season multiplicity. Each selected season retains all its May–November storm and non-storm days. Both model components used the same product weights, and counterfactual predictions were re-averaged over the resampled prediction population. Intervals are the 2.5th and 97.5th percentiles of estimable draws.

Season blocks retain shared multi-day and cross-ZIP storm shocks, as well as other within-season shared shocks; ZIP resampling accommodates dependence within ZIP across seasons. This is ZIP-and-season inference, not a bootstrap of individual storms. It assumes independence for observations sharing neither a ZIP nor a season, and the small number of seasons still limits accuracy.

Valid draws at ≥34, ≥55 and ≥64 kt: 999, 999, 999. Failed/non-estimable draws are retained in failed_draws.csv and are not silently replaced. Draw diagnostics distinguish absent phase variation from other fitting failures.

The ≥64 kt sample contains ten storms across six seasons. Limited independent hurricane information remains; 999 resamples do not add independent events. Leave-one-storm-out and leave-one-season-out influence checks for ≥34 kt are available in output/leave_one_out_34/README.md; higher-threshold checks have not been run.

| Threshold | Phase | Estimate | 95% two-way interval | Valid draws |
|---|---|---:|---:|---:|
| ≥34 kt | Cleanup | 0.33 | [-0.01, 0.58] | 999 |
| ≥34 kt | Impact | -1.50 | [-2.33, -0.84] | 999 |
| ≥34 kt | Threat | 0.29 | [-0.14, 0.63] | 999 |
| ≥55 kt | Cleanup | 0.61 | [-0.07, 1.36] | 999 |
| ≥55 kt | Impact | -2.86 | [-4.59, -1.11] | 999 |
| ≥55 kt | Threat | 0.27 | [-0.22, 0.59] | 999 |
| ≥64 kt | Cleanup | 1.04 | [0.14, 2.36] | 999 |
| ≥64 kt | Impact | -3.96 | [-7.96, -2.02] | 999 |
| ≥64 kt | Threat | 0.51 | [-0.41, 0.84] | 999 |

## Reproduce

```sh
Rscript scripts/regenerate_suppression_figures.R --cluster=zip_season --output=output/twoway_bootstrap --daily=false --mean=3.8 --bootstrap=999 --workers=8
Rscript scripts/plot_twoway_figure1.R
Rscript scripts/write_suppression_supplement.R
```

The implementation is checked against explicit copying of independently resampled ZIP and season clusters (product_weight_validation.csv). The weighted daily-data compression is also checked against the full daily-data model. New point estimates match the saved ZIP-only point estimates to numerical precision.

Historical ZIP-only intervals remain in output/suppression_corrected/Figure1_estimates.csv and comparison_with_zip_only.csv. The prior PNG is archived under figs/archive_zip_only_inference_2026-10-19/. Figure S1 remains an exploratory daily plot without sampling confidence intervals and was not refitted by this analysis.

Methodological reference: Davezies, D'Haultfoeuille and Guyonvarch, Asymptotic results under multiway clustering, https://arxiv.org/abs/1807.07925. The finite-season analysis here remains approximate.
