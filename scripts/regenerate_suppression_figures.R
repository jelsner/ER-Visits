# Confirmed suppression: absent ZIP-days are interval-observed at 0--5;
# reported counts >=6 are exact. Run from the repository root.
# Options: --bootstrap=200 --workers=4 --mean=3.8 --thresholds=34,55,64
suppressPackageStartupMessages({library(data.table);library(fixest);library(lubridate)})
setDTthreads(1); setFixest_nthreads(1)
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default) {
  x <- grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(x)) sub(paste0("^--",name,"="),"",x[1]) else default
}
B <- as.integer(option("bootstrap","200")); workers <- as.integer(option("workers","4"))
central_mean <- as.numeric(option("mean","3.8"))
thresholds <- as.integer(strsplit(option("thresholds","34,55,64"),",")[[1]])
stopifnot(B>=0,workers>=1,central_mean>=0,central_mean<=5)
out <- "output/suppression_corrected"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
note <- function(...) {cat(format(Sys.time(),"%H:%M:%S"),...,"\n");flush.console()}
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
fwrite(data.table(suppressed_lower=0,suppressed_upper=5,reported_min=6,
                 central_mean=central_mean,bootstrap_requested=B,workers=workers,
                 seed=20261018),file.path(out,"configuration.csv"))

# Reuse the previously constructed, unchanged RMW/ZIP geometry and date pairs.
pairs <- fread("output/event_time/zip_storm_alignment.csv",colClasses=c(zip="character"))
pairs[,zip:=sprintf("%05d",as.integer(zip))]
pairs[,`:=`(shared=as.IDate(shared),local=as.IDate(local))]
population <- fread("output/event_time/population_annual.csv",colClasses=c(zip="character"))
population[,zip:=sprintf("%05d",as.integer(zip))]
stopifnot(!anyDuplicated(pairs[,.(zip,SID)]),!anyDuplicated(population[,.(zip,year)]))
raw <- fread("data/Combined_ER_Data.csv")
stopifnot(!anyNA(raw$Total_Visits),all(raw$Total_Visits>=6),
          all(raw$Total_Visits==as.integer(raw$Total_Visits)))
raw <- raw[ZipCode>0 & ZipCode<1e6]
raw[,zip:=sprintf("%05d",as.integer(ZipCode))]
raw <- raw[nchar(zip)==5]
raw[,date:=as.IDate(Date)]
duplicate_rows <- sum(duplicated(raw[,.(zip,date)]))
# Retain the original pipeline's aggregation of multiple reported rows per day.
ed <- raw[,.(er_visits=sum(Total_Visits)),by=.(zip,date)]
fwrite(data.table(observed_records=nrow(ed),reported_sixes=sum(ed$er_visits==6),
                 additional_reported_rows=duplicate_rows,
                 minimum_reported=min(ed$er_visits),maximum_reported=max(ed$er_visits)),
       file.path(out,"source_audit.csv"))
rm(raw);gc()

make_calendar <- function(pp,alignment) {
  pp <- copy(pp);pp[,anchor:=get(alignment)]
  x <- pp[,.(date=as.IDate(anchor+(-8:8)),day=-8:8),by=.(zip,SID)]
  x[,priority:=fifelse(day==0,3L,fifelse(day<0,2L,1L))]
  overlaps <- x[,.(storms=uniqueN(SID)),by=.(zip,date)][storms>1]
  x[,distance:=abs(day)]
  setorderv(x,c("zip","date","priority","distance","SID"),c(1,1,-1,1,1))
  list(calendar=x[,.SD[1],by=.(zip,date)][,.(zip,date,day,SID)],overlaps=overlaps)
}

# Exact sufficient-statistic compression: covariates are constant within cells.
# The weighted binomial mean and weighted Poisson mean have the same scores
# and fixed-effect estimates as the uncompressed daily observations.
aggregate_cells <- function(dat) {
  z <- dat[,.(n=.N,ns=sum(suppressed),no=sum(!suppressed),
               sy=sum(er_visits,na.rm=TRUE),pop=first(pop)),
           by=.(zip,year_month,dow,event)]
  z[,`:=`(suppressed_fraction=ns/n,observed_mean=fifelse(no>0,sy/pmax(no,1),NA_real_),
            log_pop=log(pop))]
  z[,event:=factor(event,levels=levels(dat$event))]
  z
}
fit_pair <- function(g,multiplicity=NULL) {
  d <- copy(g)
  if(is.null(multiplicity)) d[,mult:=1] else {
    d[,mult:=unname(multiplicity[zip])];d <- d[mult>0]
  }
  d[,`:=`(wa=n*mult,wb=no*mult)]
  a <- suppressWarnings(feglm(suppressed_fraction ~ i(event,ref="None") |
      zip+year_month+dow,data=d,family=binomial("logit"),weights=~wa,notes=FALSE))
  b <- fepois(observed_mean ~ i(event,ref="None") | zip+year_month+dow,
               data=d[no>0],weights=~wb,offset=~log_pop,notes=FALSE)
  stopifnot(isTRUE(a$convStatus),isTRUE(b$convStatus))
  target <- d[,.(n=sum(n),pop=first(pop),mult=first(mult),log_pop=first(log_pop)),
              by=.(zip,year_month,dow)]
  target[,event:=factor("None",levels=levels(g$event))]
  la <- predict(a,newdata=target,type="link")
  mu <- predict(b,newdata=target,type="response")
  list(a=a,b=b,target=target,la=la,mu=mu)
}
gcomp <- function(m,levels,mask=NULL) {
  ok <- is.finite(m$la)&is.finite(m$mu)
  if(!is.null(mask)) ok <- ok & mask
  stopifnot(any(ok))
  w <- with(m$target[ok],n*pop*mult)
  avg <- function(x) sum(w*x)/sum(w)
  p0 <- plogis(m$la[ok]);mu0 <- m$mu[ok]
  baseA <- avg(p0);baseB <- avg((1-p0)*mu0)
  vals <- lapply(levels,function(lev) {
    nm <- paste0("event::",lev)
    aa <- unname(coef(m$a)[nm]);bb <- unname(coef(m$b)[nm])
    stopifnot(is.finite(aa),is.finite(bb))
    p <- plogis(m$la[ok]+aa)
    data.table(event=lev,A=avg(p-p0),C=avg((1-p)*mu0*exp(bb)-(1-p0)*mu0),
      baseA=baseA,baseC=baseB,standardization_zips=uniqueN(m$target$zip[ok]),
      standardization_zip_days=sum(m$target$n[ok]))
  })
  rbindlist(vals)
}
expand_means <- function(x) {
  rbindlist(lapply(sort(unique(c(0,2.5,3.8,5,central_mean))),function(mm) {
    z <- copy(x);z[,`:=`(mean_suppressed=mm,difference=C+mm*A,
                        baseline=baseC+mm*baseA)]
    z[,`:=`(expected=baseline+difference,percent=100*difference/baseline)];z
  }))
}
compact_model <- function(m) list(
  logit=list(coefficients=coef(m$a),fixed_effects=fixef(m$a),converged=m$a$convStatus),
  poisson=list(coefficients=coef(m$b),fixed_effects=fixef(m$b),converged=m$b$convStatus))

for(thr in thresholds) {
  note("Preparing threshold",thr)
  pp <- pairs[Wind>=thr];zz <- sort(unique(pp$zip))
  dates <- as.IDate(seq(as.Date("2005-01-01"),as.Date("2022-12-31"),by="day"))
  dates <- dates[month(dates)%in%5:11]
  base <- merge(CJ(zip=zz,date=dates),ed,by=c("zip","date"),all.x=TRUE)
  # Never replace these unknown counts with zero in the model data.
  base[,suppressed:=is.na(er_visits)]
  base[,`:=`(year=year(date),year_month=format(date,"%Y-%m"),dow=wday(date))]
  base <- merge(base,population,by=c("zip","year"),all.x=TRUE,sort=FALSE)
  stopifnot(all(is.finite(base$pop)))
  npop <- sum(base$pop<=0);base[pop<=0,pop:=1]
  stopifnot(all(base[suppressed==FALSE,er_visits]>=6),
            all(is.na(base[suppressed==TRUE,er_visits])))
  calendars <- lapply(c("shared","local"),function(a) make_calendar(pp,a))
  names(calendars) <- c("shared","local")
  phase_dat <- copy(base);phase_dat[,event:="None"]
  cc <- calendars$shared$calendar
  cc[,phase:=fifelse(day==0,"Impact",fifelse(day<0,"Threat","Cleanup"))]
  phase_dat[cc,on=.(zip,date),event:=i.phase]
  phase_dat[,event:=factor(event,levels=c("None","Threat","Impact","Cleanup"))]
  suppression_audit <- phase_dat[,.(zip_days=.N,absent_0_to_5=sum(suppressed),
    reported_sixes=sum(er_visits==6,na.rm=TRUE),suppression_fraction=mean(suppressed)),by=event]
  suppression_audit[,threshold:=thr]
  fwrite(suppression_audit,sprintf("%s/suppression_audit_%s.csv",out,thr))
  g <- aggregate_cells(phase_dat)
  note("Phase fit",thr,";",nrow(base),"daily rows compressed to",nrow(g),"cells")
  started <- proc.time()[3];m <- fit_pair(g)
  phase <- gcomp(m,c("Threat","Impact","Cleanup"))
  phase[,`:=`(threshold=thr,alignment="shared")]
  fwrite(expand_means(phase),sprintf("%s/phase_estimates_%s.csv",out,thr))
  saveRDS(compact_model(m),sprintf("%s/phase_fit_%s.rds",out,thr))
  note("Phase fit and predictions took",round(proc.time()[3]-started,1),"seconds")
  # Equivalence check for one threshold against a full daily-data fit.
  if(thr==34 && !file.exists(file.path(out,"compression_validation.csv"))) {
    note("Checking compressed model against original daily rows")
    checkA <- suppressWarnings(feglm(suppressed ~ i(event,ref="None") |
      zip+year_month+dow,data=phase_dat,family=binomial("logit"),notes=FALSE))
    checkB <- fepois(er_visits ~ i(event,ref="None") |zip+year_month+dow,
      data=phase_dat[suppressed==FALSE],offset=~log(pop),notes=FALSE)
    errA <- max(abs(coef(m$a)-coef(checkA)))
    errB <- max(abs(coef(m$b)-coef(checkB)))
    stopifnot(errA<1e-6,errB<1e-6)
    fwrite(data.table(max_logit_coefficient_error=errA,max_count_coefficient_error=errB),
           file.path(out,"compression_validation.csv"))
    rm(checkA,checkB);gc()
    note("Compression agrees with full data to",signif(max(errA,errB),3))
  }
  rm(phase_dat,m);gc()

  # Resample entire ZIP histories, as in original Figure 1. Each replicate fits
  # both components jointly on the same ZIP draw and redoes standardization.
  boot_file <- sprintf("%s/bootstrap_%s.csv",out,thr)
  existing <- if(file.exists(boot_file)) fread(boot_file) else data.table()
  completed <- if(nrow(existing)) unique(existing$replicate) else integer()
  todo <- setdiff(seq_len(B),completed)
  one_boot <- function(i) {
    tryCatch({
      set.seed(20261018+thr*1000+i)
      mult <- table(factor(sample(zz,length(zz),replace=TRUE),levels=zz))
      fit <- fit_pair(g,mult)
      v <- gcomp(fit,c("Threat","Impact","Cleanup"))
      v[,`:=`(replicate=i,threshold=thr)];v
    },error=function(e) data.table(replicate=i,error=conditionMessage(e),threshold=thr))
  }
  if(length(todo)) for(start in seq(1,length(todo),by=workers*2)) {
    ids <- todo[start:min(length(todo),start+workers*2-1)]
    started <- proc.time()[3]
    rr <- parallel::mclapply(ids,one_boot,mc.cores=workers,mc.preschedule=FALSE,
                            mc.set.seed=FALSE)
    existing <- rbindlist(c(list(existing),rr),fill=TRUE)
    fwrite(existing,boot_file)
    note("Bootstrap",thr,":",length(unique(existing$replicate)),"/",B,
         "; batch seconds",round(proc.time()[3]-started,1))
  }
  rm(g);gc()

  daily_models <- list();daily_counts <- list()
  for(align in c("shared","local")) {
    d <- copy(base);d[,event:="None"]
    cal <- calendars[[align]]$calendar
    d[cal,on=.(zip,date),event:=sprintf("d%+03d",i.day)]
    d[,event:=factor(event,levels=c("None",sprintf("d%+03d",-8:8)))]
    daily_counts[[align]] <- d[,.N,by=event]
    gg <- aggregate_cells(d);rm(d);gc()
    note("Daily fit",thr,align)
    daily_models[[align]] <- fit_pair(gg)
    rm(gg);gc()
  }
  # Both targets have the same covariates and rows, independent of event labels.
  mm <- daily_models
  keys <- c("zip","year_month","dow")
  ix <- do.call(paste,c(mm$shared$target[,..keys],sep="|"))
  iy <- do.call(paste,c(mm$local$target[,..keys],sep="|"))
  order_local <- match(ix,iy);stopifnot(!anyNA(order_local))
  valid_shared <- is.finite(mm$shared$la)&is.finite(mm$shared$mu)
  valid_local <- is.finite(mm$local$la)&is.finite(mm$local$mu)
  mask_shared <- valid_shared & valid_local[order_local]
  mask_local <- mask_shared[match(iy,ix)]
  dd <- rbindlist(lapply(c("shared","local"),function(align) {
    mask <- if(align=="shared") mask_shared else mask_local
    z <- gcomp(mm[[align]],sprintf("d%+03d",-8:8),mask)
    z[,`:=`(day=as.integer(sub("d","",event)),threshold=thr,alignment=align)]
    z <- merge(z,daily_counts[[align]],by="event",all.x=TRUE)
    setnames(z,"N","observed_zip_days");expand_means(z)
  }))
  fwrite(dd,sprintf("%s/daily_estimates_%s.csv",out,thr))
  saveRDS(lapply(daily_models,compact_model),sprintf("%s/daily_fits_%s.rds",out,thr))
  audit <- rbindlist(lapply(c("shared","local"),function(align) data.table(
    threshold=thr,alignment=align,zip_days=nrow(base),zips=uniqueN(base$zip),
    storms=uniqueN(calendars[[align]]$calendar[date%in%dates,SID]),
    overlapping_zip_days=nrow(calendars[[align]]$overlaps[date%in%dates]),
    absent_interval_days=sum(base$suppressed),nonpositive_population_guard_days=npop)))
  fwrite(audit,sprintf("%s/sample_audit_%s.csv",out,thr))
  rm(base,daily_models,mm,dd);gc()
  note("Completed threshold",thr)
}
note("Completed all requested fits. Render with scripts/plot_suppression_figures.R")
