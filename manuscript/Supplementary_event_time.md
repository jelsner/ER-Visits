# Supplementary information: corrected suppression and impact-date sensitivity

## Confirmed suppression definition and prediction assumptions

The data provider confirmed that counts of one through five visits were suppressed. Zero-count ZIP-days were also absent from the extract. An absent ZIP-day was therefore treated as an interval-observed count between zero and five; reported counts of six or more were treated as observed. Absent records were not assigned an observed count of zero, and reported sixes were not labeled as suppressed.

The logistic component modeled the probability of an absent (0–5 visit) ZIP-day, p. The count component modeled the conditional mean m among observed ZIP-days with at least six visits, using a Poisson working mean model. Combined expected visits were calculated as c p + (1-p)m. The central figures use c = 3.8 visits, David Hsu's estimate from comparison with county totals. The provider confirmed the interval 0–5; it did not establish this mean. Sensitivity analyses varied a common c across the full interval from zero to five. This assumes the suppressed-group mean is constant across ZIP-days and phases for each sensitivity scenario; these ranges are not identification bounds allowing phase-specific means to vary independently.

The count component includes ZIP, year-month, and day-of-week fixed effects and an annual ZIP population offset. The logistic component uses the same fixed effects. Population-weighted counterfactual predictions were averaged over ZIP-days with finite component predictions. As in the original pipeline, multiple reported rows within a ZIP-day were summed. The available extract supports May–November 2005–2022; the 2004 storm season is absent.

## Figure 1: corrected phase estimates

![Figure 1. Corrected phase estimates](../figs/Fig1_phase_effects_by_intensity.png)

**Figure 1. Phase-specific ED utilization around tropical cyclones with two-way ZIP-and-season inference.** Population-weighted differences in expected ED visits per ZIP-day relative to non-storm conditions, by phase and nested intensity threshold. The RMW swaths, shared storm dates, point estimates, and suppressed-group mean of 3.8 are retained. Blue bars show 95% percentile intervals from 999 independent resamples of ZIPs and the 18 calendar-year hurricane seasons. Both model components use the same product of ZIP and season multiplicities; combined predictions are re-averaged over the resampled population. Intervals use estimable draws; valid counts at ≥34, ≥55 and ≥64 kt are 999, 999, 999, respectively. Failed/non-estimable draws are recorded and not silently replaced. Orange ranges show sensitivity to a common suppressed-group mean from zero to five and are not confidence intervals. The season blocks retain shared storm shocks across ZIPs and days, including multiple storms within a year. The small number of independent hurricane seasons still limits inference.

Differences in expected visits per ZIP-day, with 95% two-way ZIP-and-season bootstrap intervals in brackets:

| Threshold | Threat | Impact | Cleanup | Non-storm baseline |
|---|---:|---:|---:|---:|
| ≥34 kt | 0.29 [-0.14, 0.63] | -1.50 [-2.33, -0.84] | 0.33 [-0.01, 0.58] | 16.44 |
| ≥55 kt | 0.27 [-0.22, 0.59] | -2.86 [-4.59, -1.11] | 0.61 [-0.07, 1.36] | 16.45 |
| ≥64 kt | 0.51 [-0.41, 0.84] | -3.96 [-7.96, -2.02] | 1.04 [0.14, 2.36] | 16.56 |

## Figure S1: daily event-time sensitivity

![Figure S1. Daily event-time alignment](../figs/FigS1_daily_event_time_alignment.png)

**Figure S1. Daily ED utilization by impact-date alignment under the confirmed suppression definition.** Daily differences in population-weighted expected ED visits per ZIP-day relative to non-storm conditions, from eight days before through eight days after the specified impact date. Panels show nested ≥34, ≥55, and ≥64 kt storm thresholds. Blue circles use a common date for all exposed ZIP areas within each storm: the date of the first retained IBTrACS track segment whose RMW buffer intersects Florida. Orange triangles use the date of closest approach to an interior representative point of each ZIP Code Tabulation Area. The vertical band marks day zero. Both alignments retain the same RMW swaths and ZIP-area intersection rule. Lines use a suppressed-group mean of 3.8; shading shows sensitivity to a common mean from zero through five. No sampling confidence intervals are shown in Figure S1. Dates use UTC.

The three phase indicators were replaced with 17 daily indicators spanning days -8 through +8. ZIP-days outside all selected storms' event windows were the reference; day -1 was not the reference. Separate models were fitted by threshold and alignment. For local timing, each ZIP area's interior representative point was projected onto retained qualifying track segments, the closest location was selected, and its timestamp was interpolated linearly between successive track fixes. This is a geographic timing proxy, not the time of maximum wind at individual patient locations or a population-weighted ZIP center.

Within each threshold, both alignment curves use the same ZIP-day prediction population. These populations contain 671, 667, 630 ZIP areas at ≥34, ≥55, and ≥64 kt, respectively. There are 39, 19, 10 contributing storms at those thresholds. Differences across threshold panels remain descriptive because their target populations differ and their storm samples are nested.

Overlapping event windows follow the original Impact > Threat > Cleanup priority, with ties within a phase resolved by the event day closest to zero and then storm identifier. At ≥34 kt, ten ZIP-days overlap under shared-date alignment and eight under local-date alignment; higher thresholds have none. Windows of storms excluded by an intensity threshold can enter that threshold's reference category, retaining the original analysis design.

The following summaries use the central suppressed-group mean. Three-day values are equal-weight averages of daily estimates over days -1 through +1, not coefficients from a refitted three-day-window model.

| Threshold | Shared date: day zero | ZIP-specific date: day zero | Shared date: three-day mean | ZIP-specific date: three-day mean |
|---|---:|---:|---:|---:|
| ≥34 kt | -1.50 | -1.07 | -0.52 | -0.46 |
| ≥55 kt | -2.86 | -2.22 | -1.08 | -0.93 |
| ≥64 kt | -3.96 | -2.69 | -1.37 | -1.02 |

## Suggested Methods wording

> Daily ZIP-code counts of one through five ED visits were suppressed by the data provider. Because zero-count ZIP-days were also absent, an absent ZIP-day represented an unknown count from zero to five. Reported counts of six or more were treated as observed. We modeled the probability of an absent record using logistic regression and modeled counts among observed records using a Poisson conditional-mean model. Combined expected visits used a suppressed-group mean of 3.8, estimated from comparison with county totals, with sensitivity analyses varying this common mean from zero to five.

## Reproducibility and remaining author checks

- Corrected fits and bootstrap: `scripts/regenerate_suppression_figures.R`.
- Figure rendering: `scripts/plot_suppression_figures.R`.
- Numerical outputs and audits: `output/suppression_corrected/`.
- Primary Figure 1 inference: `output/twoway_bootstrap/`; rendering: `scripts/plot_twoway_figure1.R`.
- Figures made before the correction are archived in `figs/archive_before_suppression_correction_2026-10-18/`.

Before submission, document the county-total calculation supporting 3.8, verify ED date conventions, reconcile main-text years and event counts, and acknowledge the limited independent hurricane seasons. Leave-one-storm-out and leave-one-season-out influence checks for ≥34 kt are available in output/leave_one_out_34/README.md; higher-threshold influence checks have not been run. Historical ZIP-only intervals are retained in output/suppression_corrected/; the primary Figure 1 uses output/twoway_bootstrap/.
