# Corrected suppression figures

Absent ZIP-days represent 0–5 visits; reported counts >=6 are observed.

Central suppressed-group mean: 3.8. Sensitivity: common mean from 0 to 5.

Run from the repository root:

```sh
Rscript scripts/regenerate_suppression_figures.R --bootstrap=200 --workers=4 --mean=3.8
Rscript scripts/plot_suppression_figures.R
Rscript scripts/write_suppression_supplement.R
```

The first command resumes saved ZIP-bootstrap replicates. To change the central mean without refitting, change central_mean in configuration.csv and rerun the final two commands; every output includes results for means 0, 2.5, 3.8, and 5. Bootstrap replicates store the linear components needed for any common mean.

Model fits use sufficient-statistic compression within identical ZIP/year-month/day-of-week/event cells. compression_validation.csv compares the compressed model with the full daily-data fit. All absent outcomes remain NA in the model data; observed sixes enter the count component.

The primary Figure 1 now uses 999 attempted two-way ZIP-and-season bootstrap draws per threshold; see output/twoway_bootstrap/ for the primary intervals, diagnostics and reproduction command. The 200-draw ZIP-only intervals in this directory are retained as historical results. Figure S1 remains a daily point-estimate and suppressed-mean sensitivity plot without confidence bands.

The RMW geometry, date definitions, population series, and nested threshold design are reused from the earlier analysis. Both corrected figures use the available May–November 2005–2022 ED record.

See manuscript/Supplementary_event_time.md for captions, methods, numerical results, and limitations.
