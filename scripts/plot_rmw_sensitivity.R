# Compare original Figure 1 with the 1.5-RMW exposure sensitivity.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out<-"output/rmw_sensitivity"
fitdir<-file.path(out,"model_1p5")
cfg<-fread(file.path(fitdir,"configuration.csv"))
stopifnot(cfg$bootstrap_mode=="zip_season",cfg$central_mean==3.8,
          cfg$end_year==2022,cfg$bootstrap_requested==999)
original<-fread("output/twoway_bootstrap/Figure1_estimates.csv")
original[,radius:="RMW"]
points<-rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/phase_estimates_%s.csv",fitdir,t))))[mean_suppressed==3.8]
draws<-rbindlist(lapply(c(34,55,64),function(t)
  fread(sprintf("%s/bootstrap_%s.csv",fitdir,t))),fill=TRUE)
if(!"error"%in%names(draws))draws[,error:=NA_character_]
draws[is.na(error)&(!is.finite(A)|!is.finite(C)),error:="non_finite_prediction"]
diagnostics<-draws[,.(attempted=uniqueN(replicate),
  valid=uniqueN(replicate[is.na(error)]),failed=uniqueN(replicate[!is.na(error)])),by=threshold]
stopifnot(nrow(diagnostics)==3,all(diagnostics$attempted==999),all(diagnostics$valid>=.95*999))
valid<-draws[is.na(error)]
stopifnot(!anyDuplicated(valid[,.(threshold,replicate,event)]),
          all(valid[,.N,by=.(threshold,replicate)]$N==3))
sm<-lapply(unique(valid$season_multiplicities),function(s)as.integer(strsplit(s,",",fixed=TRUE)[[1]]))
stopifnot(all(vapply(sm,function(s)length(s)==18 && all(s>=0) && sum(s)==18,logical(1))))
valid[,difference:=C+3.8*A]
ci<-valid[,.(ci_lo=quantile(difference,.025),ci_hi=quantile(difference,.975),
             n_boot=uniqueN(replicate)),by=.(threshold,event)]
expanded<-merge(points,ci,by=c("threshold","event"));expanded[,radius:="1.5 × RMW"]
combined<-rbindlist(list(original,expanded),fill=TRUE)
stopifnot(nrow(combined)==18,all(is.finite(combined$difference)),
          all(combined$ci_lo<combined$ci_hi))
fwrite(expanded,file.path(out,"expanded_estimates.csv"))
fwrite(combined,file.path(out,"FigureS2_estimates.csv"))
fwrite(diagnostics,file.path(out,"draw_diagnostics.csv"))
fwrite(draws[!is.na(error)],file.path(out,"failed_draws.csv"))
comparison<-merge(original[,.(threshold,event,original=difference,original_percent=percent,
                               original_ci_lo=ci_lo,original_ci_hi=ci_hi,original_support_zips=standardization_zips)],
  expanded[,.(threshold,event,expanded=difference,expanded_percent=percent,
               expanded_ci_lo=ci_lo,expanded_ci_hi=ci_hi,expanded_support_zips=standardization_zips)],
  by=c("threshold","event"))
comparison[,change:=expanded-original]
comparison[,same_sign:=sign(expanded)==sign(original)]
fwrite(comparison,file.path(out,"comparison.csv"))
sample_comparison<-rbindlist(lapply(c(34,55,64),function(t) {
  a<-fread(sprintf("output/twoway_bootstrap/sample_audit_%s.csv",t))
  b<-fread(sprintf("%s/sample_audit_%s.csv",fitdir,t))
  data.table(threshold=t,original_model_zips=a$zips,expanded_model_zips=b$zips,
    original_overlap_days=a$overlapping_zip_days,expanded_overlap_days=b$overlapping_zip_days,
    original_zip_days=a$zip_days,expanded_zip_days=b$zip_days)
}))
fwrite(sample_comparison,file.path(out,"sample_comparison.csv"))
phase_counts<-rbindlist(lapply(c(34,55,64),function(t) {
  a<-fread(sprintf("output/twoway_bootstrap/suppression_audit_%s.csv",t))
  b<-fread(sprintf("%s/suppression_audit_%s.csv",fitdir,t))
  merge(a[,.(threshold,event,original_zip_days=zip_days)],
        b[,.(threshold,event,expanded_zip_days=zip_days)],by=c("threshold","event"))
}))
fwrite(phase_counts,file.path(out,"phase_assignment_comparison.csv"))
combined[,event:=factor(event,levels=c("Threat","Impact","Cleanup"))]
combined[,radius:=factor(radius,levels=c("RMW","1.5 × RMW"))]
combined[,threshold_label:=factor(paste0("≥",threshold," kt"),levels=paste0("≥",c(34,55,64)," kt"))]
caption<-paste0("Points: suppressed-group mean 3.8 visits; bars: 95% percentile intervals from 999 two-way ZIP–season draws per radius and threshold.\n",
  "Both model components use the same resampling weights within each draw. Original track segments, storm dates and intensity classifications are held fixed.\n",
  "Wider buffers add exposed ZIP-days and can change reference assignments and prediction populations; differences between radii are descriptive.")
pd<-position_dodge(width=.35)
baseplot<-ggplot(combined,aes(event,difference,color=radius,shape=radius))+
  geom_hline(yintercept=0,color="grey65",linewidth=.4)+
  geom_errorbar(aes(ymin=ci_lo,ymax=ci_hi),position=pd,width=.12,linewidth=.7)+
  geom_point(position=pd,size=2.7)+facet_wrap(~threshold_label,nrow=1)+
  scale_color_manual(values=c("RMW"="#1864ab","1.5 × RMW"="#cf6c26"))+
  scale_shape_manual(values=c(16,17))+
  labs(x=NULL,y="Difference in expected visits per ZIP-day",color=NULL,shape=NULL,
       subtitle="Florida, May–November 2005–2022; two-way ZIP–season bootstrap",caption=caption)+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold",size=13),legend.position="top",
    plot.caption=element_text(hjust=0,size=8.5),plot.title=element_text(face="bold",size=16),
    plot.margin=margin(14,16,14,14))
p<-baseplot+labs(title="Figure S2. Sensitivity to the width of the RMW exposure corridor")
ggsave("figs/FigS2_rmw_sensitivity.png",p,width=12,height=6.3,dpi=300,device=ragg::agg_png,bg="white")
# Separate expanded-radius reproduction of Figure 1; retain the main figure.
q<-ggplot(combined[radius=="1.5 × RMW"],aes(event,difference))+
  geom_hline(yintercept=0,color="grey65",linewidth=.4)+
  geom_linerange(aes(ymin=pmin(C,C+5*A),ymax=pmax(C,C+5*A)),color="#d97730",linewidth=5,alpha=.5)+
  geom_errorbar(aes(ymin=ci_lo,ymax=ci_hi),width=.15,color="#1864ab",linewidth=.8)+
  geom_point(color="#1864ab",size=2.7)+facet_wrap(~threshold_label,nrow=1)+
  labs(title="Phase-specific ED utilization with 1.5 × RMW exposure",
    subtitle="Florida, May–November 2005–2022; original storm timing and intensity retained",
    x=NULL,y="Difference in expected visits per ZIP-day",
    caption="Points: suppressed-group mean 3.8; blue bars: 95% two-way ZIP–season bootstrap intervals (999 draws).\nOrange ranges: sensitivity to a common suppressed-group mean of 0–5; these ranges are not confidence intervals.")+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold",size=13),plot.title=element_text(face="bold",size=16),
    plot.caption=element_text(hjust=0,size=9),plot.margin=margin(14,16,14,14))
ggsave("figs/Fig1_phase_effects_rmw_1p5.png",q,width=12,height=6.1,dpi=300,device=ragg::agg_png,bg="white")
fmt<-function(x)sprintf("%.2f",x)
impact<-comparison[event=="Impact"][order(threshold)]
response<-sprintf(paste0("We assessed sensitivity to exposure width by multiplying each retained segment's RMW by 1.5 while holding the original storm inventory, intensity classifications and dates fixed. ",
  "Impact estimates were %s, %s and %s visits per ZIP-day at ≥34, ≥55 and ≥64 kt, respectively, compared with %s, %s and %s under the original RMW definition; Threat and Cleanup estimates moved closer to zero at all three thresholds. ",
  "Figure S2 presents both radii with 95%% two-way ZIP-and-season bootstrap intervals (999 draws per radius and threshold), with predictions standardized over each radius-specific sample."),
  fmt(impact$expanded[1]),fmt(impact$expanded[2]),fmt(impact$expanded[3]),
  fmt(impact$original[1]),fmt(impact$original[2]),fmt(impact$original[3]))
writeLines(response,"manuscript/Reviewer_response_RMW_sensitivity.md")
rows<-vapply(seq_len(nrow(comparison)),function(i){z<-comparison[i];sprintf(
  "| ≥%s kt | %s | %s [%s, %s] | %s [%s, %s] |",z$threshold,z$event,
  fmt(z$original),fmt(z$original_ci_lo),fmt(z$original_ci_hi),
  fmt(z$expanded),fmt(z$expanded_ci_lo),fmt(z$expanded_ci_hi))},character(1))
figure_caption<-paste0("**Figure S2. Sensitivity of phase-specific emergency department utilization estimates to the RMW exposure radius.** ",
  "Points show population-weighted differences in modeled expected ED visits per ZIP-day relative to non-storm conditions for Threat (days −8 to −1), Impact (day 0), and Cleanup (days +1 to +8). ",
  "Panels represent nested storm-intensity thresholds of ≥34, ≥55, and ≥64 kt. Blue circles use the original radius-of-maximum-wind (RMW) buffers; orange triangles use buffers with radii multiplied by 1.5. ",
  "The same original retained track segments, storm inventory, shared dates, and intensity classifications are used for both radii; ZIP-area intersections define exposure. ",
  "Bars indicate 95% percentile confidence intervals from 999 two-way bootstrap resamples of ZIP codes and the 18 May–November seasons (2005–2022) per radius and threshold, with both model components refitted using identical resampling weights within a draw. ",
  "Models include ZIP, year-month, and day-of-week fixed effects and a population offset in the count component. Absent ZIP-days represent unknown counts of 0–5 visits, reported counts ≥6 are observed, and expected visits use a suppressed-group mean of 3.8. ",
  "Predictions are standardized over each radius-specific population with finite component predictions. The horizontal line denotes zero difference. Inter-radius differences are descriptive and reflect changes in exposure assignments and prediction populations.")
methods<-c("Original segment RMW values were imputed with the original within-storm, pressure-group, within-storm sequence before Florida clipping. The exact original Florida-intersecting segments were then buffered using RMW × 1.5 × 1852 metres (RMW is recorded in nautical miles), unioned within storm and intersected with the same 2020 ZIP-area polygons. This is a controlled widening of the original corridor; it does not add storms or track segments or redefine the event date. Storm classification remains the maximum wind among the original retained segments.",
  "Each radius uses its own exposure-defined ZIP universe and phase assignments. The wider corridor can reclassify reference ZIP-days as exposed, create additional overlapping windows, and add ZIPs to the model and prediction population. The original Impact > Threat > Cleanup overlap priority and nearest-day tie rule are retained. Annual population for added ZIPs was reconstructed from the same cached ACS data, with interpolation and endpoint carry-forward. Values for original ZIPs match exactly.",
  "Both models were refitted using the same two-part conditional-mean specification as Figure 1. The expanded-radius analysis uses 999 independent ZIP-and-season product-weight draws per threshold; seeds follow the original bootstrap scheme. No paired bootstrap test of the inter-radius contrast was computed; matching draw numbers do not guarantee matching resampled seasons when the ZIP universes differ. The individual confidence intervals should not be interpreted as a formal test of the difference between radii. The finite number of seasons limits uncertainty estimation.")
exposure<-fread(file.path(out,"exposure_comparison.csv"))
exposure_rows<-vapply(seq_len(nrow(exposure)),function(i){z<-exposure[i];sprintf(
  "| ≥%s kt | %s | %s | %s | %s | %s |",z$threshold,z$storms,z$original_zips,z$expanded_zips,z$original_pairs,z$expanded_pairs)},character(1))
writeLines(c("# Figure S2: RMW-radius sensitivity","",figure_caption,"",
  "## Methods","",methods,"","## Results","",response,"",
  "Differences in expected visits per ZIP-day; brackets give 95% two-way ZIP-and-season bootstrap intervals.","",
  "| Threshold | Phase | RMW | 1.5 × RMW |","|---|---|---:|---:|",rows,"",
  "## Exposure inventory (2005–2022)","",
  "| Threshold | Storms | Original ZIPs | Expanded ZIPs | Original ZIP–storm pairs | Expanded ZIP–storm pairs |",
  "|---|---:|---:|---:|---:|---:|",exposure_rows),"manuscript/Supplementary_RMW_sensitivity.md")
writeLines(c("# RMW-radius sensitivity","",figure_caption,"",methods,"",
  "| Threshold | Phase | RMW | 1.5 × RMW |","|---|---|---:|---:|",rows,"",
  "## Reproduce","","```sh","Rscript scripts/prepare_rmw_sensitivity.R",
  "Rscript scripts/regenerate_suppression_figures.R --cluster=zip_season --output=output/rmw_sensitivity/model_1p5 --pairs=output/rmw_sensitivity/zip_storm_pairs_rmw_1p5.csv --population=output/rmw_sensitivity/population_annual.csv --daily=false --mean=3.8 --bootstrap=999 --workers=8",
  "Rscript scripts/plot_rmw_sensitivity.R","```","",
  "Original RMW estimates and intervals are read from output/twoway_bootstrap/Figure1_estimates.csv. The original Figure 1 and Figure S1 are retained. Figure S2 is in figs/FigS2_rmw_sensitivity.png; the standalone expanded-radius figure is figs/Fig1_phase_effects_rmw_1p5.png.","",
  "geometry_validation.txt checks exact reproduction of the original exposure assignments, nesting and population equality. draw_diagnostics.csv and failed_draws.csv retain attempted/valid counts and fitting failures."),file.path(out,"README.md"))
protected<-fread(file.path(out,"protected_fingerprints.csv"))
protected[,current_md5:=unname(tools::md5sum(file))]
stopifnot(all(protected$md5==protected$current_md5))
fwrite(protected,file.path(out,"protected_validation.csv"))
writeLines(c("Protected original exposure/population inputs and Figures 1 and S1 are byte-identical to their pre-analysis versions.",
  paste("Valid expanded-radius bootstrap draws:",paste(diagnostics$valid,collapse=", ")),
  "Every valid draw has a complete phase triplet and season multiplicities summing to 18.",
  "Expanded point estimates and intervals are finite; original estimates and intervals are reused without modification."),
  file.path(out,"validation_summary.txt"))
print(comparison)
