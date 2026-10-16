# Daily event-time alignment diagnostic

Run from the repository root:

```sh
Rscript scripts/fetch_event_time_population.R
Rscript scripts/daily_event_time.R
```

To redraw the figure without refitting: `Rscript scripts/plot_daily_event_time.R`.

The population step downloads public ACS data and caches it in this directory.
The analysis reads existing 2020 ZCTA/state shapefiles from the local tigris
cache and uses the repository's IBTrACS archive and ED extract. It does not alter
the original manuscript, analysis notebook, or existing figures.

## What the comparison estimates

The diagnostic retains the original RMW-buffered track geometry, missing-RMW
imputation, ZIP intersection rule, nested storm thresholds, logistic model of
the suppression indicator, Poisson model for counts above six, population
offset, and ZIP/year-month/day-of-week fixed effects. Daily event labels from
-8 through +8 replace the three phase indicators. "None" is the reference:
days outside all selected storms' event windows, rather than day -1.

The shared-date curve uses `first(Date)` among qualifying track segments whose
RMW buffers intersect Florida, as in the posted code. The local-date curve uses
closest approach to an interior representative point of each ZCTA. Closest
approach is calculated by projecting that point onto each retained storm track
segment and linearly interpolating its timestamp between successive fixes.
This is a geographic timing proxy, not a local wind maximum or a population-
weighted ZIP location. Dates use UTC to match the existing date construction;
the ED date convention should be verified before manuscript use.

Each curve is obtained from a separately fitted daily model. Predictions are
combined using `6*p + (1-p)*m`, preserving the original suppressed-group mean
assumption. Predictions are population-weighted over the same rows with finite
predictions from both components of both alignment models within each wind
threshold. Logistic fixed effects can drop ZIPs with no suppression variation;
the resulting common prediction support is reported in `daily_estimates.csv`.
Target populations differ across threshold panels, so comparisons across panels
remain descriptive rather than formal tests of intensity differences.

Overlaps follow the original Impact > Threat > Cleanup priority. Within-phase
ties choose the event day nearest zero, then storm SID. Counts of overlapping
ZIP-days appear in `sample_audit.csv`. Lower-intensity storms excluded by a
threshold can enter that threshold's "None" reference, matching the original
pipeline; correcting this is a separate robustness analysis.

## Limits and interpretation

The available ED data do not include the 2004 hurricane season. Models use
May-November 2005-2022. The storm inventory still contains the original 44
swaths, but only events with windows in the model period contribute exposure.
The sample audit reports contributing storm counts separately.

The combined curves show point estimates only, with no confidence bands. They
are a diagnostic of timing and magnitude, not a significance test. The component
coefficient file includes two-way ZIP/year clustered intervals (18 year clusters)
for the logistic and Poisson components; these are not confidence intervals for
the combined visit estimates. Joint uncertainty accounting for storm dependence
is needed before inferential manuscript use.

The missing-row-to-zero assignment and suppressed-group mean of six are retained
to isolate the alignment comparison. Those assumptions still need provider
confirmation and sensitivity analysis. Population estimates available in the
download begin in 2011; earlier model years carry the earliest available value
within each ZIP, using the posted interpolation/extrapolation rule.

`timing_summary.csv` reports day zero, the equal-weight average for days -1 to
+1, pre/post eight-day averages, and the minimum daily estimate. These averages
are summaries of the daily model, not coefficients from a refitted phase model.

## Outputs

- `../../figs/daily_event_time_alignment.png`: display figure.
- `daily_estimates.csv`: combined daily contrasts, baselines, percentages, support.
- `timing_summary.csv`: timing and magnitude summaries.
- `zip_storm_alignment.csv`: paired shared/local dates and shifts.
- `storm_inventory.csv`: original swath inventory.
- `sample_audit.csv`: model-period event counts, overlaps, zero assignments.
- `component_coefficients.csv`: component estimates and diagnostic intervals.
- `models_*.rds`: compact coefficients, fixed effects, convergence and sample summaries.
- `population_annual.csv`: annual interpolated/extrapolated population inputs.
- `session_info.txt`: R/package versions.
