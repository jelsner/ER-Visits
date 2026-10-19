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

The corrected models use the available May–November 2005–2022 ED record. Figure
1 retains the original 200-resample ZIP-bootstrap uncertainty design; Figure S1
shows daily point estimates and suppressed-mean sensitivity, without sampling
confidence intervals. RMW exposure geometry and timing definitions are retained.
