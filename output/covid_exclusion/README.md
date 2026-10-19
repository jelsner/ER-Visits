# Sensitivity analysis excluding 2020–2022

The corrected phase models were refitted using May–November 2005–2019, excluding every observation from 2020, 2021, and 2022. The RMW exposure definition, nested intensity thresholds, phase timing, fixed effects, population offset/weighting, and suppressed-group mean of 3.8 were retained. Reported counts >=6 are observed; absent ZIP-days represent 0–5 visits. Each threshold used 200 ZIP-cluster bootstrap resamples with both components refitted.

Estimates are differences in expected visits per ZIP-day relative to non-storm conditions; brackets show 95% ZIP-bootstrap intervals.

| Threshold | Phase | Full period (2005–2022) | Excluding 2020–2022 (2005–2019) |
|---|---|---:|---:|
| ≥34 kt | Threat | 0.29 [0.24, 0.35] | 0.08 [0.01, 0.14] |
| ≥34 kt | Impact | -1.50 [-1.65, -1.35] | -1.52 [-1.67, -1.37] |
| ≥34 kt | Cleanup | 0.33 [0.28, 0.40] | 0.29 [0.22, 0.35] |
| ≥55 kt | Threat | 0.27 [0.16, 0.40] | 0.08 [-0.04, 0.20] |
| ≥55 kt | Impact | -2.86 [-3.12, -2.60] | -3.09 [-3.44, -2.73] |
| ≥55 kt | Cleanup | 0.61 [0.47, 0.76] | 0.80 [0.61, 1.02] |
| ≥64 kt | Threat | 0.51 [0.32, 0.68] | 0.19 [-0.02, 0.40] |
| ≥64 kt | Impact | -3.96 [-4.28, -3.64] | -5.70 [-6.17, -5.03] |
| ≥64 kt | Cleanup | 1.04 [0.86, 1.21] | 2.10 [1.78, 2.40] |

Predictions are standardized within each fitted period, on ZIP-days with finite predictions from both model components. Consequently, these comparisons include changes in the observed period, storm sample, population weights, and prediction support; they are not formal tests of COVID interactions or a causal attribution of estimate changes to COVID.

The existing ZIP-cluster inference design is retained; this analysis does not resolve the reviewers' separate concern about dependence across ZIPs affected by the same storm.

## Reproduce

```sh
Rscript scripts/regenerate_suppression_figures.R --end-year=2019 --output=output/covid_exclusion --daily=false --mean=3.8 --bootstrap=200 --workers=4
Rscript scripts/summarize_covid_exclusion.R
```
