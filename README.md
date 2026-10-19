# ER-Visits

## Corrected Figures 1 and S1

The provider confirmed that absent ZIP-days represent unknown counts from zero
through five; reported counts of six or more are observed. The corrected figures
use a central suppressed-group mean of 3.8 visits (David Hsu's estimate from a
county-total comparison), with sensitivity to a common mean from zero to five.

Run from the repository root:

```sh
Rscript scripts/regenerate_suppression_figures.R --bootstrap=200 --workers=4 --mean=3.8
Rscript scripts/plot_suppression_figures.R
Rscript scripts/write_suppression_supplement.R
```

The fitting script resumes saved ZIP-bootstrap replicates. The rendering and
supplement-writing steps can be repeated without refitting.

- Current figures: `figs/Fig1_phase_effects_by_intensity.png` and
  `figs/FigS1_daily_event_time_alignment.png`.
- Corrected estimates, bootstrap draws, audits and methods:
  `output/suppression_corrected/`.
- Draft supplementary captions and results:
  `manuscript/Supplementary_event_time.md`.

The original submission notebooks, CSVs directly under `output/`, and
`output/event_time/` are retained as historical analyses and use the superseded
suppression interpretation. Earlier
figure versions and event-time scripts are archived in directories named
`archive_before_suppression_correction_2026-10-18` under `figs/` and `scripts/`.

The corrected models use the available May–November 2005–2022 ED record.
Figure S1 shows daily point estimates and suppressed-mean sensitivity, without
sampling confidence intervals. RMW exposure geometry and timing definitions
are retained.

## Primary Figure 1 inference: ZIP and season

Figure 1 now uses 999 two-way ZIP-and-season bootstrap draws per threshold.
ZIPs and the 18 May–November seasons are independently sampled with replacement;
both model components use the same product weights. Point estimates are
unchanged. The 200-draw ZIP-only intervals are retained as historical results.

```sh
Rscript scripts/regenerate_suppression_figures.R --cluster=zip_season --output=output/twoway_bootstrap --daily=false --mean=3.8 --bootstrap=999 --workers=8
Rscript scripts/plot_twoway_figure1.R
Rscript scripts/write_suppression_supplement.R
```

Primary intervals, comparison with ZIP-only intervals, draw diagnostics and
validation are in `output/twoway_bootstrap/`. The previous Figure 1 is archived
in `figs/archive_zip_only_inference_2026-10-19/`. Existing plotting entry points
retain the selected two-way inference when redrawing Figure 1.

Intervals remain approximate with only 18 seasons and ten qualifying hurricanes
across six seasons. For ≥34 kt, leave-one-storm-out (39 storms) and
leave-one-season-out (18 seasons) checks are now complete; see
`output/leave_one_out_34/README.md` and `figs/leave_one_out_34.png`.
Reproduce with `Rscript scripts/leave_one_out_34.R`. These are influence
checks without new confidence intervals; higher-threshold checks were not run.

## Exclusion of COVID years

The phase models were refitted excluding 2020–2022, retaining 2005–2019, with
the same corrected suppression treatment and central mean of 3.8. Results,
200-resample ZIP-bootstrap intervals and the full-period comparison are in
`output/covid_exclusion/README.md` and `output/covid_exclusion/comparison.csv`.
The two-sentence reviewer response is saved in
`manuscript/Reviewer_response_COVID_exclusion.md`. The existing full-period
figures and estimates are retained.

The COVID exclusion intervals still use ZIP-only resampling; they have not been
updated to the two-way design used in the primary Figure 1.
