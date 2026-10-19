# Leave-one-out influence checks: ≥34 kt

Both model components were refitted after each of 39 storm omissions and 18 season omissions (2005–2022). The ≥34 kt RMW exposure rule, shared-date phase calendar, fixed effects, population offset, reported-count cutoff ≥6 and suppressed-group mean 3.8 were retained. The ZIP universe remains the full ≥34 kt universe.

For each omitted storm, every ZIP-day in its exposed ZIPs from day −8 through +8 was removed from estimation and prediction, including overlapping windows. The remaining original phase assignments were retained. Omitted storm days were never recoded as non-storm days. Non-window observations from those ZIPs remain. A season omission removes the entire May–November period across all ZIPs, including non-storm days, even in seasons without a qualifying storm.

Population-weighted predictions are standardized over retained ZIP-days with finite predictions from both components. Consequently, estimates reflect changes in the fitted associations and the retained standardization population. These are descriptive influence checks, not new confidence intervals, significance tests, or independent replications. Figure 1 retains the full-sample two-way bootstrap intervals. Additional c=0, 2.5, and 5 estimates are saved in all_estimates.csv.

| Omission | Phase | Full estimate | Leave-one-out range | Largest influence | Absolute change |
|---|---|---:|---:|---|---:|
| Leave one storm out | Threat | 0.29 | 0.15 to 0.34 | ALEX 2022 | 0.15 |
| Leave one storm out | Impact | -1.50 | -1.79 to -1.20 | ALEX 2022 | 0.29 |
| Leave one storm out | Cleanup | 0.33 | 0.25 to 0.37 | WILMA 2005 | 0.08 |
| Leave one season out | Threat | 0.29 | 0.04 to 0.33 | 2022 | 0.25 |
| Leave one season out | Impact | -1.50 | -1.68 to -1.17 | 2017 | 0.33 |
| Leave one season out | Cleanup | 0.33 | 0.23 to 0.37 | 2005 | 0.10 |

Validation: all 58 paired fits converged; full-model contrast agreement with Figure 1: maximum absolute error 2.44e-15. Removed rows by original phase are audited in removed_days_by_phase.csv.

Reproduce: `Rscript scripts/leave_one_out_34.R`.

Plot: ../../figs/leave_one_out_34.png. Full numerical output: central_estimates.csv and all_estimates.csv.
