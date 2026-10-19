# Run from the repository root: Rscript scripts/daily_event_time.R
# Population download, if needed: Rscript scripts/fetch_event_time_population.R
suppressPackageStartupMessages({
  library(sf); library(data.table); library(dplyr); library(lubridate)
  library(fixest); library(ggplot2)
})
setFixest_nthreads(4)
dir.create("output/event_time", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
log_note <- function(...) { cat(..., "\n"); flush.console() }
log_note("Preparing RMW swaths and ZIP-specific timing")

# Reproduce the posted RMW construction, including imputation before FL clipping.
tracks <- st_read("data/IBTrACS.NA.list.v04r00.lines.shp", quiet = TRUE) |>
  st_transform(32616) |>
  filter(year >= 2004, year <= 2022, USA_WIND >= 34) |>
  group_by(SID) |>
  mutate(USA_RMW = ifelse(is.na(USA_RMW), mean(USA_RMW, na.rm=TRUE), USA_RMW)) |>
  group_by(USA_PRES) |>
  mutate(USA_RMW = ifelse(is.na(USA_RMW), mean(USA_RMW, na.rm=TRUE), USA_RMW)) |>
  group_by(SID) |>
  mutate(USA_RMW = ifelse(is.na(USA_RMW), mean(USA_RMW, na.rm=TRUE), USA_RMW)) |>
  ungroup()
fl <- USAboundaries::us_states(states="FL", resolution="low") |> st_transform(32616)
buffers <- st_buffer(tracks, dist=tracks$USA_RMW * 1852)
keep <- lengths(st_intersects(buffers, fl)) > 0
tracks <- tracks[keep, ]; buffers <- buffers[keep, ]
buffers$Date <- as.Date(buffers$ISO_TIME)
swaths <- buffers |> group_by(SID) |>
  summarize(Date0=first(Date), NAME=first(NAME), Wind=max(USA_WIND),
            geometry=st_union(geometry), .groups="drop")
stopifnot(nrow(swaths) == 44)

# Read the same 2020 cartographic ZCTA boundaries from the user's existing cache.
geo_cache <- path.expand("~/Library/Caches/tigris")
zips <- st_read(file.path(geo_cache,"cb_2020_us_zcta520_500k.shp"),quiet=TRUE)
states <- st_read(file.path(geo_cache,"cb_2020_us_state_500k.shp"),quiet=TRUE)
zips <- zips[lengths(st_intersects(zips, states[states$STUSPS=="FL",])) > 0,]
zips$zip <- zips$ZCTA5CE20
zips <- st_transform(zips,32616)
# Keep the geographic intersection used in the manuscript's ZIP assignment.
hits <- st_intersects(st_transform(zips,4326), st_transform(swaths,4326))
pairs <- rbindlist(lapply(seq_along(hits), function(i) {
  j <- hits[[i]]
  if (!length(j)) return(NULL)
  data.table(zip=zips$zip[i], SID=swaths$SID[j], NAME=swaths$NAME[j],
             Wind=swaths$Wind[j], shared=swaths$Date0[j], zi=i)
}))

# Local alignment: closest approach to a point guaranteed to lie inside each ZCTA.
# Project each point onto each retained line segment. The timestamp at the next
# vertex is normally three hours later in this archive; derive it from the next
# ISO_TIME within the original storm track, rather than assume a fixed interval.
all_tracks <- st_read("data/IBTrACS.NA.list.v04r00.lines.shp", quiet=TRUE)
times <- as.data.table(st_drop_geometry(all_tracks))[,.(SID,ISO_TIME)]
times[, time := as.POSIXct(ISO_TIME,tz="UTC")]
setorder(times,SID,time)
times[, end_time := shift(time,type="lead"),by=SID]
times <- unique(times,by=c("SID","ISO_TIME"))
track_dt <- as.data.table(st_drop_geometry(tracks))
track_dt <- merge(track_dt,times[,.(SID,ISO_TIME,time,end_time)],
                  by=c("SID","ISO_TIME"),sort=FALSE)
# Preserve exact row order for the geometry-to-time lookup.
lookup <- paste(track_dt$SID,track_dt$ISO_TIME)
track_dt <- track_dt[match(paste(tracks$SID,tracks$ISO_TIME),lookup)]
points <- st_coordinates(suppressWarnings(st_point_on_surface(zips)))
pairs[, local := as.Date(NA)]
pairs[, local_time_utc := as.POSIXct(NA,tz="UTC")]
pairs[, closest_distance_km := NA_real_]
for (r in seq_len(nrow(pairs))) {
  ix <- which(tracks$SID == pairs$SID[r])
  p <- points[pairs$zi[r],1:2]
  candidates <- lapply(ix, function(j) {
    xy <- st_coordinates(tracks[j,])[,1:2,drop=FALSE]
    a <- xy[1,]; b <- xy[nrow(xy),]; v <- b-a
    f <- if(sum(v*v)>0) max(0,min(1,sum((p-a)*v)/sum(v*v))) else 0
    dist <- sqrt(sum((p-(a+f*v))^2))
    t0 <- as.numeric(track_dt$time[j]); t1 <- as.numeric(track_dt$end_time[j])
    if(!is.finite(t1)) t1 <- t0
    c(dist=dist,time=t0+f*(t1-t0))
  })
  mat <- do.call(rbind,candidates); best <- which.min(mat[,"dist"])
  tm <- as.POSIXct(mat[best,"time"],origin="1970-01-01",tz="UTC")
  set(pairs,r,"local_time_utc",tm)
  set(pairs,r,"local",as.Date(tm,tz="UTC"))
  set(pairs,r,"closest_distance_km",mat[best,"dist"]/1000)
}
pairs[, shift_days := as.integer(local-shared)]
fwrite(pairs[,!"zi"],"output/event_time/zip_storm_alignment.csv")
fwrite(as.data.table(st_drop_geometry(swaths)),"output/event_time/storm_inventory.csv")
log_note("Swaths:",nrow(swaths),"; ZIP-storm pairs:",nrow(pairs))

# Preserve source count preparation and missing-row-to-zero assumption.
raw <- fread("data/Combined_ER_Data.csv")
raw <- raw[ZipCode > 0 & ZipCode < 1e6]
raw[,zip := sprintf("%05d",as.integer(ZipCode))]
raw <- raw[nchar(zip)==5]
raw[,date := as.IDate(Date)]
ed <- raw[,.(er_visits=sum(Total_Visits,na.rm=TRUE)),by=.(zip,date)]
rm(raw, all_tracks); gc()

# ACS interpolation/extrapolation matches the manuscript's posted approach.
pop_files <- list.files("output/event_time",pattern="^acs_[0-9]+[.]rds$",full.names=TRUE)
stopifnot(length(pop_files)>0)
annual <- rbindlist(lapply(pop_files,function(f) {
  z <- as.data.table(readRDS(f))
  z[,.(zip=GEOID,year=as.integer(sub(".*acs_([0-9]+)[.]rds","\\1",f)),pop=estimate)]
}))[zip %in% pairs$zip]
population <- annual[,.(year=2005:2022,
  pop=approx(year,pop,xout=2005:2022,method="linear",rule=2)$y),by=zip]
stopifnot(!anyDuplicated(population[,.(zip,year)]))
fwrite(population,"output/event_time/population_annual.csv")

# Resolve overlaps using original phase priority, then nearest relative day and
# SID to make ties within a phase deterministic. Save overlaps for transparency.
make_calendar <- function(pp,alignment) {
  pp <- copy(pp); pp[,anchor := get(alignment)]
  expanded <- pp[,.(date=as.IDate(anchor+(-8:8)),event_day=-8:8),
                 by=.(zip,SID,Wind)]
  expanded[,priority := fifelse(event_day==0,3L,fifelse(event_day<0,2L,1L))]
  overlaps <- expanded[,.(n_storms=uniqueN(SID)),by=.(zip,date)][n_storms>1]
  expanded[,distance := abs(event_day)]
  setorderv(expanded,c("zip","date","priority","distance","SID"),c(1,1,-1,1,1))
  resolved <- expanded[, .SD[1],by=.(zip,date)]
  list(calendar=resolved[,.(zip,date,event_day,SID)],overlaps=overlaps)
}

results <- list(); components <- list(); audits <- list(); summaries <- list()
for(thr in c(34L,55L,64L)) {
  pp <- pairs[Wind>=thr]
  zz <- sort(unique(pp$zip))
  dates <- as.IDate(seq(as.Date("2005-01-01"),as.Date("2022-12-31"),by="day"))
  dates <- dates[month(dates)%in%5:11]
  base <- merge(CJ(zip=zz,date=dates),ed,by=c("zip","date"),all.x=TRUE)
  base[,imputed_zero := is.na(er_visits)]
  base[is.na(er_visits),er_visits:=0L]
  base[,`:=`(year=year(date),year_month=format(date,"%Y-%m"),dow=wday(date),
             is_le6=as.integer(er_visits==6L))]
  base <- merge(base,population,by=c("zip","year"),all.x=TRUE,sort=FALSE)
  stopifnot(all(is.finite(base$pop)))
  # Match the original guard for unpopulated ZCTAs; report these explicitly.
  n_pop_guard <- sum(base$pop<=0)
  base[pop<=0,pop:=1]
  base[,log_pop:=log(pop)]
  models <- list(); datasets <- list(); common <- rep(TRUE,nrow(base))
  for(alignment in c("shared","local")) {
    cal <- make_calendar(pp,alignment)
    dat <- copy(base)
    dat[,event := "None"]
    dat[cal$calendar,on=.(zip,date),event:=sprintf("d%+03d",i.event_day)]
    dat[,event := factor(event,levels=c("None",sprintf("d%+03d",-8:8)))]
    log_note("Fitting",thr,"kt",alignment,"on",nrow(dat),"ZIP-days")
    a <- feglm(is_le6 ~ i(event,ref="None") | zip+year_month+dow,
               data=dat,family=binomial("logit"),notes=FALSE)
    b <- fepois(er_visits ~ i(event,ref="None") | zip+year_month+dow,
                data=dat[er_visits>6],offset=~log_pop,notes=FALSE)
    stopifnot(isTRUE(a$convStatus),isTRUE(b$convStatus))
    nd <- copy(dat); nd[,event:=factor("None",levels=levels(dat$event))]
    la <- predict(a,newdata=nd,type="link")
    mb <- predict(b,newdata=nd,type="response")
    common <- common & is.finite(la) & is.finite(mb)
    models[[alignment]] <- list(a=a,b=b,la=la,mb=mb)
    datasets[[alignment]] <- dat
    audits[[paste(thr,alignment)]] <- data.table(threshold=thr,alignment=alignment,
      zip_days=nrow(dat),zips=uniqueN(dat$zip),
      storms=uniqueN(cal$calendar[date%in%dates,SID]),
      overlapping_zip_days=nrow(cal$overlaps[date%in%dates]),
      source_missing_zero_days=sum(dat$imputed_zero),
      nonpositive_population_guard_days=n_pop_guard)
  }
  stopifnot(sum(common)>0)
  for(alignment in names(models)) {
    m <- models[[alignment]]; dat <- datasets[[alignment]]
    w <- dat$pop[common]
    average <- function(x) weighted.mean(x,w)
    baseline <- average(6*plogis(m$la[common])+(1-plogis(m$la[common]))*m$mb[common])
    res <- rbindlist(lapply(-8:8,function(k) {
      nm <- paste0("event::",sprintf("d%+03d",k))
      ca <- coef(m$a)[nm]; cb <- coef(m$b)[nm]
      stopifnot(is.finite(ca),is.finite(cb))
      p <- plogis(m$la[common]+ca)
      ey <- average(6*p+(1-p)*m$mb[common]*exp(cb))
      data.table(threshold=thr,alignment=alignment,day=k,
                 baseline=baseline,expected=ey,difference=ey-baseline,
                 percent=100*(ey-baseline)/baseline,
                 observed_zip_days=sum(dat$event==sprintf("d%+03d",k)),
                 standardization_zip_days=sum(common),
                 standardization_zips=uniqueN(dat$zip[common]))
    }))
    results[[paste(thr,alignment)]] <- res
    # Component intervals are diagnostics only: two-way ZIP/year clustering.
    # Year clusters retain statewide cross-ZIP and cross-day dependence, with
    # only 18 years; these are not joint intervals for the combined curve.
    for(part in c("a","b")) {
      ct <- as.data.table(coeftable(m[[part]],vcov=~zip+year),keep.rownames="term")
      ci <- confint(m[[part]],vcov=~zip+year)
      ct[,`:=`(threshold=thr,alignment=alignment,component=part,
                ci_low=ci[,1],ci_high=ci[,2])]
      components[[paste(thr,alignment,part)]] <- ct
    }
    # Summary is an equal-weight average of day estimates, not a refit.
    summaries[[paste(thr,alignment)]] <- res[,.(
      day_zero=difference[day==0],
      days_minus1_to_plus1=mean(difference[day%in%-1:1]),
      pre_8_days=mean(difference[day<0]),post_8_days=mean(difference[day>0]),
      trough_day=day[which.min(difference)],trough_difference=min(difference)
    ),by=.(threshold,alignment)]
  }
  # Keep reproducible fit summaries without multi-gigabyte observation arrays.
  fit_summary <- lapply(models,function(m) list(
    logit=list(coefficients=coef(m$a),fixed_effects=fixef(m$a),
               converged=m$a$convStatus,nobs=nobs(m$a)),
    poisson=list(coefficients=coef(m$b),fixed_effects=fixef(m$b),
                 converged=m$b$convStatus,nobs=nobs(m$b))))
  saveRDS(list(models=fit_summary,standardization_rows=sum(common)),
          sprintf("output/event_time/models_%s.rds",thr))
  rm(base,datasets,models,dat); gc()
}
res <- rbindlist(results); audit <- rbindlist(audits); summary <- rbindlist(summaries)
fwrite(res,"output/event_time/daily_estimates.csv")
fwrite(audit,"output/event_time/sample_audit.csv")
fwrite(summary,"output/event_time/timing_summary.csv")
fwrite(rbindlist(components,fill=TRUE),"output/event_time/component_coefficients.csv")
source("scripts/plot_daily_event_time.R")
