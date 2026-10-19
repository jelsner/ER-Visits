# Compare the corrected full-period phase models with the 2005--2019 refits.
suppressPackageStartupMessages(library(data.table))
out <- "output/covid_exclusion"
cfg <- fread(file.path(out,"configuration.csv"))
stopifnot(cfg$end_year==2019,cfg$central_mean==3.8,cfg$include_daily==FALSE)
full <- fread("output/suppression_corrected/Figure1_estimates.csv")
stopifnot(all(full$mean_suppressed==3.8))
pre <- rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/phase_estimates_%s.csv",out,t))))[mean_suppressed==3.8]
boot <- rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/bootstrap_%s.csv",out,t))),fill=TRUE)
stopifnot(nrow(pre)==9,nrow(boot)==1800,
          all(boot[,uniqueN(replicate),by=threshold]$V1==200),
          !anyDuplicated(boot[,.(threshold,replicate,event)]),
          all(is.finite(boot$A)),all(is.finite(boot$C)))
boot[,difference:=C+3.8*A]
ci <- boot[,.(ci_lo=quantile(difference,.025),ci_hi=quantile(difference,.975),
               n_boot=uniqueN(replicate)),by=.(threshold,event)]
pre <- merge(pre,ci,by=c("threshold","event"))
fwrite(pre,file.path(out,"phase_estimates_with_ci.csv"))
comparison <- merge(
  full[,.(threshold,event,full_estimate=difference,full_ci_lo=ci_lo,full_ci_hi=ci_hi,
           full_percent=percent,full_baseline=baseline,full_support_zips=standardization_zips)],
  pre[,.(threshold,event,excluded_estimate=difference,excluded_ci_lo=ci_lo,
          excluded_ci_hi=ci_hi,excluded_percent=percent,excluded_baseline=baseline,
          excluded_support_zips=standardization_zips)],by=c("threshold","event"))
comparison[,change_in_estimate:=excluded_estimate-full_estimate]
comparison[,event:=factor(event,levels=c("Threat","Impact","Cleanup"))]
setorder(comparison,threshold,event)
fwrite(comparison,file.path(out,"comparison.csv"))
audits <- rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/sample_audit_%s.csv",out,t))))
fwrite(audits,file.path(out,"sample_audit.csv"))
for(t in c(34,55,64)) {
  m <- readRDS(sprintf("%s/phase_fit_%s.rds",out,t))
  for(component in m) {
    stopifnot(component$converged)
    years <- as.integer(substr(names(component$fixed_effects$year_month),1,4))
    stopifnot(length(years)>0,min(years)==2005,max(years)==2019,
              !any(years%in%2020:2022))
  }
}
fmt <- function(x) sprintf("%.2f",x)
rows <- apply(comparison,1,function(r) sprintf("| ≥%s kt | %s | %s [%s, %s] | %s [%s, %s] |",
  r[["threshold"]],r[["event"]],fmt(as.numeric(r[["full_estimate"]])),
  fmt(as.numeric(r[["full_ci_lo"]])),fmt(as.numeric(r[["full_ci_hi"]])),
  fmt(as.numeric(r[["excluded_estimate"]])),fmt(as.numeric(r[["excluded_ci_lo"]])),
  fmt(as.numeric(r[["excluded_ci_hi"]]))))
writeLines(c("# Sensitivity analysis excluding 2020–2022","",
"The corrected phase models were refitted using May–November 2005–2019, excluding every observation from 2020, 2021, and 2022. The RMW exposure definition, nested intensity thresholds, phase timing, fixed effects, population offset/weighting, and suppressed-group mean of 3.8 were retained. Reported counts >=6 are observed; absent ZIP-days represent 0–5 visits. Each threshold used 200 ZIP-cluster bootstrap resamples with both components refitted.","",
"Estimates are differences in expected visits per ZIP-day relative to non-storm conditions; brackets show 95% ZIP-bootstrap intervals.","",
"| Threshold | Phase | Full period (2005–2022) | Excluding 2020–2022 (2005–2019) |",
"|---|---|---:|---:|",rows,"",
"Predictions are standardized within each fitted period, on ZIP-days with finite predictions from both model components. Consequently, these comparisons include changes in the observed period, storm sample, population weights, and prediction support; they are not formal tests of COVID interactions or a causal attribution of estimate changes to COVID.","",
"The existing ZIP-cluster inference design is retained; this analysis does not resolve the reviewers' separate concern about dependence across ZIPs affected by the same storm.","",
"## Reproduce","","```sh",
"Rscript scripts/regenerate_suppression_figures.R --end-year=2019 --output=output/covid_exclusion --daily=false --mean=3.8 --bootstrap=200 --workers=4",
"Rscript scripts/summarize_covid_exclusion.R","```"),file.path(out,"README.md"))
writeLines(c("Six point-model components converged; all fixed-effect month labels fall in 2005–2019.",
  "600 successful ZIP-bootstrap replicates; no duplicate replicate/threshold/phase rows.",
  "Nine finite phase contrasts, with finite bootstrap intervals.",
  "Suppression definition, central mean, exposure and model structure retained."),
  file.path(out,"validation_summary.txt"))
print(comparison);print(audits)
