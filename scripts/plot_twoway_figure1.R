# Publish Figure 1 using the ZIP x season (pigeonhole) bootstrap.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- "output/twoway_bootstrap"
cfg <- fread(file.path(out,"configuration.csv"))
stopifnot(cfg$bootstrap_mode=="zip_season",cfg$end_year==2022,cfg$central_mean==3.8)
mm <- cfg$central_mean[1];B<-cfg$bootstrap_requested[1]
phase <- rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/phase_estimates_%s.csv",out,t))))[mean_suppressed==mm]
boot <- rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/bootstrap_%s.csv",out,t))),fill=TRUE)
if(!"error"%in%names(boot)) boot[,error:=NA_character_]
boot[is.na(error)&(!is.finite(A)|!is.finite(C)),error:="non_finite_prediction"]
diagnostics <- boot[,.(attempted=uniqueN(replicate),
  valid=uniqueN(replicate[is.na(error)&is.finite(A)&is.finite(C)]),
  missing_phase=uniqueN(replicate[!is.na(error)&
    grepl("non_estimable_missing_phase",error)]),
  other_failures=uniqueN(replicate[!is.na(error)&
    !grepl("non_estimable_missing_phase",error)])),by=threshold]
setorder(diagnostics,threshold)
stopifnot(nrow(diagnostics)==3,all(diagnostics$attempted==B),
          all(diagnostics$valid>=.95*B))
valid <- boot[is.na(error)&is.finite(A)&is.finite(C)]
stopifnot(!anyDuplicated(valid[,.(threshold,replicate,event)]),
          all(valid[,.N,by=.(threshold,replicate)]$N==3))
fwrite(boot[!is.na(error)],file.path(out,"failed_draws.csv"))
fwrite(diagnostics,file.path(out,"draw_diagnostics.csv"))
valid[,difference:=C+mm*A]
ci <- valid[,.(ci_lo=quantile(difference,.025),ci_hi=quantile(difference,.975),
                n_boot=uniqueN(replicate)),by=.(threshold,event)]
phase <- merge(phase,ci,by=c("threshold","event"))
phase[,`:=`(mean_lo=pmin(C,C+5*A),mean_hi=pmax(C,C+5*A),
             n_attempted=B,inference="ZIP x hurricane season")]
old <- fread("output/suppression_corrected/Figure1_estimates.csv")
comparison <- merge(old[,.(threshold,event,estimate_zip=difference,
                           ci_lo_zip=ci_lo,ci_hi_zip=ci_hi)],
                     phase[,.(threshold,event,estimate_two_way=difference,
                               ci_lo_two_way=ci_lo,ci_hi_two_way=ci_hi)],
                     by=c("threshold","event"))
point_err <- max(abs(comparison$estimate_zip-comparison$estimate_two_way))
stopifnot(nrow(phase)==9,point_err<1e-8,all(is.finite(phase$ci_lo)),
          all(is.finite(phase$ci_hi)),all(phase$ci_lo<phase$ci_hi))
comparison[,width_ratio:=(ci_hi_two_way-ci_lo_two_way)/(ci_hi_zip-ci_lo_zip)]
fwrite(phase,file.path(out,"Figure1_estimates.csv"))
fwrite(comparison,file.path(out,"comparison_with_zip_only.csv"))
phase[,event:=factor(event,levels=c("Threat","Impact","Cleanup"))]
phase[,threshold_label:=factor(paste0("\u2265",threshold," kt"),
                 levels=paste0("\u2265",c(34,55,64)," kt"))]
draw_text <- if(all(diagnostics$valid==B))
  sprintf("%s draws per threshold",B) else
  sprintf("%s draws/threshold; valid: %s",B,paste(diagnostics$valid,collapse=", "))
caption <- sprintf(paste0(
  "Blue bars: 95%% two-way ZIP\u2013season bootstrap intervals (%s). Points use a suppressed-group mean of 3.8.\n",
  "ZIPs and the 18 May\u2013November seasons were independently resampled; both model components used the same resampling weights.\n",
  "Orange ranges: sensitivity of point estimates to a common suppressed-group mean of 0\u20135; these are not confidence intervals."),draw_text)
p <- ggplot(phase,aes(event,difference))+
  geom_hline(yintercept=0,color="#8b949e",linewidth=.4)+
  geom_linerange(aes(ymin=mean_lo,ymax=mean_hi),color="#d97730",linewidth=5,alpha=.5)+
  geom_errorbar(aes(ymin=ci_lo,ymax=ci_hi),color="#1864ab",width=.15,linewidth=.8)+
  geom_point(color="#1864ab",size=2.8)+facet_wrap(~threshold_label,nrow=1)+
  scale_y_continuous(breaks=scales::breaks_width(2))+
  labs(title="Figure 1. Phase-specific ED utilization around tropical cyclones",
       subtitle="Two-way ZIP\u2013season inference; RMW storm-core exposure, 2005\u20132022",
       x=NULL,y="Difference in expected visits per ZIP-day",caption=caption)+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold",size=13),
    plot.title=element_text(face="bold",size=17),
    plot.caption=element_text(hjust=0,size=8.5),plot.margin=margin(15,18,15,15))
archive <- "figs/archive_zip_only_inference_2026-10-19"
dir.create(archive,showWarnings=FALSE)
dest <- file.path(archive,"Fig1_phase_effects_by_intensity.png")
if(file.exists("figs/Fig1_phase_effects_by_intensity.png")&&!file.exists(dest))
  stopifnot(file.copy("figs/Fig1_phase_effects_by_intensity.png",dest))
ggsave("figs/Fig1_phase_effects_by_intensity.png",p,width=11.5,height=6.1,
       dpi=300,device=ragg::agg_png,bg="white")
fmt <- function(x) sprintf("%.2f",x)
rows <- vapply(seq_len(nrow(phase)),function(i) {
  z<-phase[i];sprintf("| ≥%s kt | %s | %s | [%s, %s] | %s |",z$threshold,z$event,
                      fmt(z$difference),fmt(z$ci_lo),fmt(z$ci_hi),z$n_boot)
},character(1))
writeLines(c("# Two-way ZIP-and-season bootstrap for Figure 1","",
"The main point estimates, corrected 0–5 suppression definition, suppressed-group mean of 3.8, RMW swaths, phase calendar, population adjustment and prediction target are retained. Only the resampling design changes.","",
sprintf("Within each threshold, %s bootstrap draws independently resampled all ZIPs and the 18 calendar-year hurricane seasons with replacement. A cell's frequency weight was multiplied by its ZIP multiplicity and season multiplicity. Each selected season retains all its May–November storm and non-storm days. Both model components used the same product weights, and counterfactual predictions were re-averaged over the resampled prediction population. Intervals are the 2.5th and 97.5th percentiles of estimable draws.",B),"",
"Season blocks retain shared multi-day and cross-ZIP storm shocks, as well as other within-season shared shocks; ZIP resampling accommodates dependence within ZIP across seasons. This is ZIP-and-season inference, not a bootstrap of individual storms. It assumes independence for observations sharing neither a ZIP nor a season, and the small number of seasons still limits accuracy.","",
sprintf("Valid draws at ≥34, ≥55 and ≥64 kt: %s. Failed/non-estimable draws are retained in failed_draws.csv and are not silently replaced. Draw diagnostics distinguish absent phase variation from other fitting failures.",paste(diagnostics$valid,collapse=", ")),"",
"The ≥64 kt sample contains ten storms across six seasons. Limited independent hurricane information remains; 999 resamples do not add independent events. Leave-one-storm-out and leave-one-season-out influence checks for ≥34 kt are available in output/leave_one_out_34/README.md; higher-threshold checks have not been run.","",
"| Threshold | Phase | Estimate | 95% two-way interval | Valid draws |",
"|---|---|---:|---:|---:|",rows,"",
"## Reproduce","","```sh",
"Rscript scripts/regenerate_suppression_figures.R --cluster=zip_season --output=output/twoway_bootstrap --daily=false --mean=3.8 --bootstrap=999 --workers=8",
"Rscript scripts/plot_twoway_figure1.R",
"Rscript scripts/write_suppression_supplement.R","```","",
"The implementation is checked against explicit copying of independently resampled ZIP and season clusters (product_weight_validation.csv). The weighted daily-data compression is also checked against the full daily-data model. New point estimates match the saved ZIP-only point estimates to numerical precision.","",
"Historical ZIP-only intervals remain in output/suppression_corrected/Figure1_estimates.csv and comparison_with_zip_only.csv. The prior PNG is archived under figs/archive_zip_only_inference_2026-10-19/. Figure S1 remains an exploratory daily plot without sampling confidence intervals and was not refitted by this analysis.","",
"Methodological reference: Davezies, D'Haultfoeuille and Guyonvarch, Asymptotic results under multiway clustering, https://arxiv.org/abs/1807.07925. The finite-season analysis here remains approximate."),
  file.path(out,"README.md"))
writeLines(c(sprintf("Maximum point-estimate change: %.3g",point_err),
             "All bootstrap draws and failures are retained; confidence intervals use valid complete phase triplets.",
             paste("Valid draws by threshold:",paste(diagnostics$valid,collapse=", "))),
           file.path(out,"validation_summary.txt"))
print(phase[,.(threshold,event,difference,ci_lo,ci_hi,n_boot)])
print(diagnostics)
