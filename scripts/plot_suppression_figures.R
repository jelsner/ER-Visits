# Render corrected Figures 1 and S1 from saved fits; no model refitting.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- "output/suppression_corrected"
cfg <- fread(file.path(out,"configuration.csv")); mm <- cfg$central_mean[1]
read_set <- function(prefix) {
  files <- sprintf("%s/%s_%s.csv",out,prefix,c(34,55,64))
  stopifnot(all(file.exists(files)));rbindlist(lapply(files,fread),fill=TRUE)
}
phase_all <- read_set("phase_estimates")
daily_all <- read_set("daily_estimates")
boots <- read_set("bootstrap")
if("error"%in%names(boots)) stopifnot(all(is.na(boots$error)))
stopifnot(all(boots[,uniqueN(replicate),by=threshold]$V1>=200))
phase <- phase_all[mean_suppressed==phase_all$mean_suppressed[1]]
phase[,`:=`(mean_suppressed=mm,difference=C+mm*A,baseline=baseC+mm*baseA)]
phase[,`:=`(expected=baseline+difference,percent=100*difference/baseline)]
phase[,`:=`(mean_lo=pmin(C,C+5*A),mean_hi=pmax(C,C+5*A))]
boots[,difference:=C+mm*A]
cis <- boots[,.(ci_lo=quantile(difference,.025),ci_hi=quantile(difference,.975),
                n_boot=uniqueN(replicate)),by=.(threshold,event)]
phase <- merge(phase,cis,by=c("threshold","event"))
fwrite(phase,file.path(out,"Figure1_estimates.csv"))
fwrite(phase_all,file.path(out,"phase_suppressed_mean_sensitivity.csv"))
fwrite(daily_all,file.path(out,"daily_suppressed_mean_sensitivity.csv"))
phase[,event:=factor(event,levels=c("Threat","Impact","Cleanup"))]
phase[,threshold_label:=factor(paste0("\u2265",threshold," kt"),
                               levels=paste0("\u2265",c(34,55,64)," kt"))]
main_theme <- theme_minimal(base_size=12)+theme(
  panel.grid.minor=element_blank(),strip.text=element_text(face="bold",size=13),
  plot.title=element_text(face="bold",size=17),
  plot.caption=element_text(hjust=0,size=9),plot.margin=margin(15,18,15,15))
caption1 <- sprintf(paste0(
  "Points: suppressed-group mean = %s; blue bars: 95%% ZIP-cluster bootstrap intervals (200 resamples).\n",
  "Orange ranges: sensitivity of point estimates to a suppressed-group mean of 0\u20135; these are not confidence intervals.\n",
  "Absent ZIP-days represent 0\u20135 visits; reported counts \u22656 are observed. Population-weighted contrasts versus non-storm conditions."),
  format(mm,trim=TRUE))
p1 <- ggplot(phase,aes(event,difference))+
  geom_hline(yintercept=0,color="#8b949e",linewidth=.4)+
  geom_linerange(aes(ymin=mean_lo,ymax=mean_hi),color="#d97730",linewidth=5,alpha=.5)+
  geom_errorbar(aes(ymin=ci_lo,ymax=ci_hi),color="#1864ab",width=.15,linewidth=.8)+
  geom_point(color="#1864ab",size=2.8)+facet_wrap(~threshold_label,nrow=1)+
  labs(title="Figure 1. Phase-specific ED utilization around tropical cyclones",
       subtitle="Confirmed suppression definition; RMW storm-core exposure, 2005\u20132022",
       x=NULL,y="Difference in expected visits per ZIP-day",caption=caption1)+main_theme
daily <- daily_all[mean_suppressed==daily_all$mean_suppressed[1]]
daily[,`:=`(mean_suppressed=mm,difference=C+mm*A,baseline=baseC+mm*baseA)]
daily[,`:=`(expected=baseline+difference,percent=100*difference/baseline)]
daily[,`:=`(mean_lo=pmin(C,C+5*A),mean_hi=pmax(C,C+5*A))]
fwrite(daily,file.path(out,"FigureS1_estimates.csv"))
daily[,threshold_label:=factor(paste0("\u2265",threshold," kt"),
                               levels=paste0("\u2265",c(34,55,64)," kt"))]
daily[,alignment_label:=factor(alignment,levels=c("shared","local"),
                 labels=c("Shared storm date","ZIP-specific closest approach"))]
pS <- ggplot(daily,aes(day,difference,color=alignment_label,fill=alignment_label,
                       shape=alignment_label))+
  annotate("rect",xmin=-.5,xmax=.5,ymin=-Inf,ymax=Inf,fill="#f2e7d4",alpha=.7)+
  geom_hline(yintercept=0,color="#8b949e",linewidth=.4)+
  geom_vline(xintercept=0,color="#999999",linewidth=.4,linetype="dotted")+
  geom_ribbon(aes(ymin=mean_lo,ymax=mean_hi),alpha=.15,color=NA)+
  geom_line(linewidth=.8)+geom_point(size=2)+facet_wrap(~threshold_label,nrow=1)+
  scale_color_manual(values=c("#1864ab","#c05621"))+
  scale_fill_manual(values=c("#1864ab","#c05621"))+
  scale_shape_manual(values=c(16,17))+
  scale_x_continuous(breaks=seq(-8,8,2))+
  labs(title="Figure S1. Daily ED utilization by impact-date alignment",
       subtitle="Confirmed suppression definition; same RMW storm-core swaths, 2005\u20132022",
       x="Days relative to the specified impact date",y="Difference in expected visits per ZIP-day",
       color=NULL,fill=NULL,shape=NULL,
       caption=sprintf(paste0(
         "Lines: suppressed-group mean = %s; shading: sensitivity to a mean of 0\u20135. No confidence intervals are shown.\n",
         "Reference: days outside all qualifying \u00b18-day windows. Same prediction population across alignments within each panel.\n",
         "Absent ZIP-days represent 0\u20135 visits; reported counts \u22656 are observed. Population-weighted; storm dates use UTC."),
         format(mm,trim=TRUE)))+main_theme+theme(legend.position="bottom")
dir.create("figs",showWarnings=FALSE)
archive <- "figs/archive_before_suppression_correction_2026-10-18"
dir.create(archive,showWarnings=FALSE)
names_to_keep <- c("Fig1_phase_effects_by_intensity.png","FigS1_daily_event_time_alignment.png",
                   "daily_event_time_alignment.png")
for(nm in names_to_keep) {
  src <- file.path("figs",nm);dest <- file.path(archive,nm)
  if(file.exists(src)&&!file.exists(dest)) stopifnot(file.copy(src,dest))
}
if(file.exists("output/twoway_bootstrap/Figure1_estimates.csv")) {
  # Preserve the selected primary inference when rerendering existing outputs.
  source("scripts/plot_twoway_figure1.R",local=new.env(parent=environment()))
} else {
  ggsave("figs/Fig1_phase_effects_by_intensity.png",p1,width=11.5,height=5.8,
         dpi=300,device=ragg::agg_png,bg="white")
}
ggsave("figs/FigS1_daily_event_time_alignment.png",pS,width=12.5,height=6.1,
       dpi=300,device=ragg::agg_png,bg="white")
ggsave("figs/daily_event_time_alignment.png",pS,width=12.5,height=6.1,
       dpi=300,device=ragg::agg_png,bg="white")
summary <- daily[,.(day_zero=difference[day==0],
  three_day_mean=mean(difference[day%in%-1:1]),
  pre_mean=mean(difference[day<0]),post_mean=mean(difference[day>0]),
  trough_day=day[which.min(difference)],trough_difference=min(difference)),
  by=.(threshold,alignment,mean_suppressed)]
fwrite(summary,file.path(out,"timing_summary.csv"))
print(phase[,.(threshold,event,difference,percent,ci_lo,ci_hi,mean_lo,mean_hi)])
print(summary)
