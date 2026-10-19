# Influence analysis, >=34 kt only. Run from the repository root.
suppressPackageStartupMessages({library(data.table);library(fixest);library(lubridate);library(ggplot2)})
setDTthreads(1);setFixest_nthreads(1)
out <- "output/leave_one_out_34"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
central_mean <- 3.8
# Reuse the exact fitting, phase-calendar and standardization implementation.
expr <- parse("scripts/regenerate_suppression_figures.R")
for(e in expr) if(is.call(e)&&identical(e[[1]],as.name("<-")) &&
  is.symbol(e[[2]])&&as.character(e[[2]])%in%
  c("make_calendar","aggregate_cells","fit_pair","gcomp","expand_means")) eval(e)
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
inputs <- c("data/Combined_ER_Data.csv","output/event_time/zip_storm_alignment.csv",
            "output/event_time/population_annual.csv","scripts/regenerate_suppression_figures.R")
fwrite(data.table(file=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,"input_fingerprints.csv"))
pairs <- fread(inputs[2],colClasses=c(zip="character"))
pairs[,zip:=sprintf("%05d",as.integer(zip))]
pairs[,`:=`(shared=as.IDate(shared),local=as.IDate(local))]
pp <- pairs[Wind>=34];zz <- sort(unique(pp$zip))
population <- fread(inputs[3],colClasses=c(zip="character"))
population[,zip:=sprintf("%05d",as.integer(zip))]
raw <- fread(inputs[1]);stopifnot(!anyNA(raw$Total_Visits),all(raw$Total_Visits>=6))
raw <- raw[ZipCode>0 & ZipCode<1e6]
raw[,zip:=sprintf("%05d",as.integer(ZipCode))];raw<-raw[nchar(zip)==5]
raw[,date:=as.IDate(Date)]
ed <- raw[,.(er_visits=sum(Total_Visits)),by=.(zip,date)]
rm(raw);gc()
dates <- as.IDate(seq(as.Date("2005-01-01"),as.Date("2022-12-31"),by="day"))
dates <- dates[month(dates)%in%5:11]
dat <- merge(CJ(zip=zz,date=dates),ed,by=c("zip","date"),all.x=TRUE)
dat[,suppressed:=is.na(er_visits)]
dat[,`:=`(year=year(date),year_month=format(date,"%Y-%m"),dow=wday(date))]
dat <- merge(dat,population,by=c("zip","year"),all.x=TRUE,sort=FALSE)
stopifnot(all(is.finite(dat$pop)));dat[pop<=0,pop:=1]
cal <- make_calendar(pp,"shared")$calendar
cal[,phase:=fifelse(day==0,"Impact",fifelse(day<0,"Threat","Cleanup"))]
dat[,event:="None"];dat[cal,on=.(zip,date),event:=i.phase]
dat[,event:=factor(event,levels=c("None","Threat","Impact","Cleanup"))]
# Include every ZIP-day in the omitted storm's window, including the few
# overlapping days whose assigned phase belongs to a different storm.
windows <- pp[,.(date=as.IDate(shared+(-8:8))),by=.(zip,SID)]
windows <- windows[date%in%dates]
storms <- unique(pp[SID%in%windows$SID,.(SID,NAME,shared)])
stopifnot(!anyDuplicated(storms$SID),nrow(storms)==39)
storms[,season:=year(shared)];setorder(storms,shared,SID)
fwrite(storms,file.path(out,"storm_inventory.csv"))
setkey(dat,zip,date)
jobs <- rbindlist(list(data.table(design="Full sample",omitted="None",label="Full sample",season=NA_integer_),
  storms[,.(design="Leave one storm out",omitted=SID,label=paste(NAME,season),season)],
  data.table(design="Leave one season out",omitted=as.character(2005:2022),
             label=as.character(2005:2022),season=2005:2022)))
one <- function(i) {
  job<-jobs[i]
  tryCatch({
    if(job$design=="Leave one storm out") {
      remove <- unique(windows[SID==job$omitted,.(zip,date)])
      d <- dat[!remove,on=.(zip,date)]
      stopifnot(nrow(d)==nrow(dat)-nrow(remove),!any(d[remove,on=.(zip,date),which=TRUE,nomatch=0L]>0L))
    } else if(job$design=="Leave one season out") {
      d<-dat[year!=job$season];stopifnot(!job$season%in%d$year)
    } else d<-dat
    removed <- dat[!d,on=.(zip,date)]
    audit <- removed[,.(removed_zip_days=.N),by=event]
    audit[,`:=`(design=job$design,omitted=job$omitted)]
    m<-fit_pair(aggregate_cells(d))
    z<-expand_means(gcomp(m,c("Threat","Impact","Cleanup")))
    z[,`:=`(design=job$design,omitted=job$omitted,label=job$label,omitted_season=job$season,
             retained_zip_days=nrow(d),removed_zip_days=nrow(dat)-nrow(d),
             logit_converged=m$a$convStatus,poisson_converged=m$b$convStatus)]
    list(estimates=z,audit=audit,error=NULL)
  },error=function(e)list(estimates=NULL,audit=NULL,error=data.table(design=job$design,
    omitted=job$omitted,error=conditionMessage(e))))
}
res<-list()
for(start in seq(1,nrow(jobs),by=4)) {
  ids<-start:min(start+3,nrow(jobs))
  res<-c(res,parallel::mclapply(ids,one,mc.cores=4,mc.preschedule=FALSE))
  cat(format(Sys.time(),"%H:%M:%S"),"Completed",length(res),"of",nrow(jobs),"fits\n");flush.console()
  saveRDS(res,file.path(out,"checkpoint.rds"))
}
errors<-rbindlist(lapply(res,`[[`,"error"),fill=TRUE)
if(nrow(errors)) {fwrite(errors,file.path(out,"failed_fits.csv"));stop("See failed_fits.csv")}
est<-rbindlist(lapply(res,`[[`,"estimates"))
audit<-rbindlist(lapply(res,`[[`,"audit"))
fwrite(audit,file.path(out,"removed_days_by_phase.csv"))
full<-est[design=="Full sample" & mean_suppressed==3.8]
primary<-fread("output/twoway_bootstrap/Figure1_estimates.csv")[threshold==34]
err<-max(abs(full$difference-primary$difference[match(full$event,primary$event)]))
stopifnot(err<1e-8,nrow(est[mean_suppressed==3.8])==58*3,
          all(est$logit_converged),all(est$poisson_converged))
# All means match their corresponding full-data comparison, not only 3.8.
reference<-est[design=="Full sample",.(event,mean_suppressed,ref=difference)]
est[reference,on=.(event,mean_suppressed),`:=`(full_difference=i.ref,change_from_full=difference-i.ref)]
fwrite(est,file.path(out,"all_estimates.csv"))
central<-est[mean_suppressed==3.8];fwrite(central,file.path(out,"central_estimates.csv"))
summary<-central[design!="Full sample",.(minimum=min(difference),maximum=max(difference),
  minimum_omission=label[which.min(difference)],maximum_omission=label[which.max(difference)],
  full_estimate=first(full_difference),largest_absolute_change=max(abs(change_from_full)),
  most_influential=label[which.max(abs(change_from_full))],
  all_same_sign=all(sign(difference)==sign(full_difference)),fits=.N),by=.(design,event)]
fwrite(summary,file.path(out,"summary.csv"))
central[,event:=factor(event,levels=c("Threat","Impact","Cleanup"))]
p<-ggplot(central[design!="Full sample"],aes(x=difference,y=reorder(label,omitted_season)))+
  geom_vline(xintercept=0,color="grey70",linewidth=.3)+
  geom_vline(data=central[design=="Full sample",.(event,difference)],
             aes(xintercept=difference),color="#d97730",linetype="dashed",linewidth=.6)+
  geom_point(color="#1864ab",size=1.5)+
  facet_grid(design~event,scales="free",space="free_y")+
  labs(title="Influence of individual storms and seasons: ≥34 kt",
    subtitle="Population-weighted phase contrasts; corrected suppression model, mean 3.8 visits",
    x="Difference in expected visits per ZIP-day",y=NULL,
    caption="Each dot is a refitted model after the labeled omission; dashed orange lines show full-sample estimates.\nStorm omissions remove exposed ZIP-days from day −8 to +8; season omissions remove all May–November days.\nThese are influence checks, without new confidence intervals. Each fit uses its retained prediction population.")+
  theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
    axis.text.y=element_text(size=7),strip.text=element_text(face="bold"),
    plot.caption=element_text(hjust=0,size=9),plot.title=element_text(face="bold"))
ggsave("figs/leave_one_out_34.png",p,width=12,height=15,dpi=220,device=ragg::agg_png,bg="white")
fmt<-function(x)sprintf("%.2f",x)
rows<-vapply(seq_len(nrow(summary)),function(i){z<-summary[i];sprintf(
  "| %s | %s | %s | %s to %s | %s | %s |",z$design,z$event,fmt(z$full_estimate),
  fmt(z$minimum),fmt(z$maximum),z$most_influential,fmt(z$largest_absolute_change))},character(1))
writeLines(c("# Leave-one-out influence checks: ≥34 kt","",
"Both model components were refitted after each of 39 storm omissions and 18 season omissions (2005–2022). The ≥34 kt RMW exposure rule, shared-date phase calendar, fixed effects, population offset, reported-count cutoff ≥6 and suppressed-group mean 3.8 were retained. The ZIP universe remains the full ≥34 kt universe.","",
"For each omitted storm, every ZIP-day in its exposed ZIPs from day −8 through +8 was removed from estimation and prediction, including overlapping windows. The remaining original phase assignments were retained. Omitted storm days were never recoded as non-storm days. Non-window observations from those ZIPs remain. A season omission removes the entire May–November period across all ZIPs, including non-storm days, even in seasons without a qualifying storm.","",
"Population-weighted predictions are standardized over retained ZIP-days with finite predictions from both components. Consequently, estimates reflect changes in the fitted associations and the retained standardization population. These are descriptive influence checks, not new confidence intervals, significance tests, or independent replications. Figure 1 retains the full-sample two-way bootstrap intervals. Additional c=0, 2.5, and 5 estimates are saved in all_estimates.csv.","",
"| Omission | Phase | Full estimate | Leave-one-out range | Largest influence | Absolute change |",
"|---|---|---:|---:|---|---:|",rows,"",
sprintf("Validation: all 58 paired fits converged; full-model contrast agreement with Figure 1: maximum absolute error %.3g. Removed rows by original phase are audited in removed_days_by_phase.csv.",err),"",
"Reproduce: `Rscript scripts/leave_one_out_34.R`.","",
"Plot: ../../figs/leave_one_out_34.png. Full numerical output: central_estimates.csv and all_estimates.csv."),file.path(out,"README.md"))
writeLines(c(sprintf("Full-model maximum difference: %.6g",err),
  "58 paired fits converged (full sample + 39 storm omissions + 18 season omissions).",
  "Storm removal includes every window ZIP-day, including overlaps; no exposure days recoded as reference.",
  "Season removal includes every May–November ZIP-day in the omitted year."),file.path(out,"validation_summary.txt"))
print(summary)
