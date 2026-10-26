# RMW-radius sensitivity

**Figure S2. Sensitivity of phase-specific emergency department utilization estimates to the RMW exposure radius.** Points show population-weighted differences in modeled expected ED visits per ZIP-day relative to non-storm conditions for Threat (days −8 to −1), Impact (day 0), and Cleanup (days +1 to +8). Panels represent nested storm-intensity thresholds of ≥34, ≥55, and ≥64 kt. Blue circles use the original radius-of-maximum-wind (RMW) buffers; orange triangles use buffers with radii multiplied by 1.5. The same original retained track segments, storm inventory, shared dates, and intensity classifications are used for both radii; ZIP-area intersections define exposure. Bars indicate 95% percentile confidence intervals from 999 two-way bootstrap resamples of ZIP codes and the 18 May–November seasons (2005–2022) per radius and threshold, with both model components refitted using identical resampling weights within a draw. Models include ZIP, year-month, and day-of-week fixed effects and a population offset in the count component. Absent ZIP-days represent unknown counts of 0–5 visits, reported counts ≥6 are observed, and expected visits use a suppressed-group mean of 3.8. Predictions are standardized over each radius-specific population with finite component predictions. The horizontal line denotes zero difference. Inter-radius differences are descriptive and reflect changes in exposure assignments and prediction populations.

Original segment RMW values were imputed with the original within-storm, pressure-group, within-storm sequence before Florida clipping. The exact original Florida-intersecting segments were then buffered using RMW × 1.5 × 1852 metres (RMW is recorded in nautical miles), unioned within storm and intersected with the same 2020 ZIP-area polygons. This is a controlled widening of the original corridor; it does not add storms or track segments or redefine the event date. Storm classification remains the maximum wind among the original retained segments.
Each radius uses its own exposure-defined ZIP universe and phase assignments. The wider corridor can reclassify reference ZIP-days as exposed, create additional overlapping windows, and add ZIPs to the model and prediction population. The original Impact > Threat > Cleanup overlap priority and nearest-day tie rule are retained. Annual population for added ZIPs was reconstructed from the same cached ACS data, with interpolation and endpoint carry-forward. Values for original ZIPs match exactly.
Both models were refitted using the same two-part conditional-mean specification as Figure 1. The expanded-radius analysis uses 999 independent ZIP-and-season product-weight draws per threshold; seeds follow the original bootstrap scheme. No paired bootstrap test of the inter-radius contrast was computed; matching draw numbers do not guarantee matching resampled seasons when the ZIP universes differ. The individual confidence intervals should not be interpreted as a formal test of the difference between radii. The finite number of seasons limits uncertainty estimation.

| Threshold | Phase | RMW | 1.5 × RMW |
|---|---|---:|---:|
| ≥34 kt | Cleanup | 0.33 [-0.01, 0.58] | 0.08 [-0.25, 0.38] |
| ≥34 kt | Impact | -1.50 [-2.33, -0.84] | -1.65 [-2.32, -1.00] |
| ≥34 kt | Threat | 0.29 [-0.14, 0.63] | 0.00 [-0.42, 0.36] |
| ≥55 kt | Cleanup | 0.61 [-0.07, 1.36] | 0.19 [-0.39, 0.99] |
| ≥55 kt | Impact | -2.86 [-4.59, -1.11] | -2.91 [-4.63, -1.52] |
| ≥55 kt | Threat | 0.27 [-0.22, 0.59] | -0.17 [-0.83, 0.34] |
| ≥64 kt | Cleanup | 1.04 [0.14, 2.36] | 0.50 [-0.71, 1.94] |
| ≥64 kt | Impact | -3.96 [-7.96, -2.02] | -3.95 [-7.40, -2.47] |
| ≥64 kt | Threat | 0.51 [-0.41, 0.84] | 0.02 [-1.29, 0.53] |

## Reproduce

```sh
Rscript scripts/prepare_rmw_sensitivity.R
Rscript scripts/regenerate_suppression_figures.R --cluster=zip_season --output=output/rmw_sensitivity/model_1p5 --pairs=output/rmw_sensitivity/zip_storm_pairs_rmw_1p5.csv --population=output/rmw_sensitivity/population_annual.csv --daily=false --mean=3.8 --bootstrap=999 --workers=8
Rscript scripts/plot_rmw_sensitivity.R
```

Original RMW estimates and intervals are read from output/twoway_bootstrap/Figure1_estimates.csv. The original Figure 1 and Figure S1 are retained. Figure S2 is in figs/FigS2_rmw_sensitivity.png; the standalone expanded-radius figure is figs/Fig1_phase_effects_rmw_1p5.png.

geometry_validation.txt checks exact reproduction of the original exposure assignments, nesting and population equality. draw_diagnostics.csv and failed_draws.csv retain attempted/valid counts and fitting failures.
