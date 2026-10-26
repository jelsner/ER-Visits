# Widen each original retained track-segment buffer by 50%, keeping the
# original storm set, intensity classification and shared event dates fixed.
suppressPackageStartupMessages({library(sf);library(data.table);library(dplyr);library(lubridate)})
setDTthreads(1)
out<-"output/rmw_sensitivity"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
note<-function(...) {cat(format(Sys.time(),"%H:%M:%S"),...,"\n");flush.console()}
protected<-c("output/event_time/zip_storm_alignment.csv","output/event_time/population_annual.csv",
             "output/twoway_bootstrap/Figure1_estimates.csv","figs/Fig1_phase_effects_by_intensity.png",
             "figs/FigS1_daily_event_time_alignment.png")
fwrite(data.table(file=protected,md5=unname(tools::md5sum(protected))),file.path(out,"protected_fingerprints.csv"))
tracks_file<-"data/IBTrACS.NA.list.v04r00.lines.shp"
note("Rebuilding original segment buffers")
tracks<-st_read(tracks_file,quiet=TRUE) |> st_transform(32616) |>
  filter(year>=2004,year<=2022,USA_WIND>=34) |> group_by(SID) |>
  mutate(USA_RMW=ifelse(is.na(USA_RMW),mean(USA_RMW,na.rm=TRUE),USA_RMW)) |>
  group_by(USA_PRES) |>
  mutate(USA_RMW=ifelse(is.na(USA_RMW),mean(USA_RMW,na.rm=TRUE),USA_RMW)) |>
  group_by(SID) |>
  mutate(USA_RMW=ifelse(is.na(USA_RMW),mean(USA_RMW,na.rm=TRUE),USA_RMW)) |> ungroup()
stopifnot(all(is.finite(tracks$USA_RMW)),all(tracks$USA_RMW>0))
fl<-USAboundaries::us_states(states="FL",resolution="low") |> st_transform(32616)
buffers<-st_buffer(tracks,dist=tracks$USA_RMW*1852)
keep<-lengths(st_intersects(buffers,fl))>0
tracks<-tracks[keep,];buffers<-buffers[keep,]
buffers$Date<-as.Date(buffers$ISO_TIME)
swaths<-buffers |> group_by(SID) |> summarize(Date0=first(Date),NAME=first(NAME),
  Wind=max(USA_WIND),geometry=st_union(geometry),.groups="drop")
stopifnot(nrow(swaths)==44)
expanded<-st_buffer(tracks,dist=tracks$USA_RMW*1852*1.5)
expanded<-expanded |> group_by(SID) |> summarize(geometry=st_union(geometry),.groups="drop")
meta<-as.data.table(st_drop_geometry(swaths))
expanded<-merge(expanded,as.data.frame(meta),by="SID",sort=FALSE)
geo_cache<-path.expand("~/Library/Caches/tigris")
zips<-st_read(file.path(geo_cache,"cb_2020_us_zcta520_500k.shp"),quiet=TRUE)
states<-st_read(file.path(geo_cache,"cb_2020_us_state_500k.shp"),quiet=TRUE)
zips<-zips[lengths(st_intersects(zips,states[states$STUSPS=="FL",]))>0,]
zips$zip<-zips$ZCTA5CE20
assign_pairs<-function(sw) {
  hits<-st_intersects(st_transform(zips,4326),st_transform(sw,4326))
  rbindlist(lapply(seq_along(hits),function(i) {
    j<-hits[[i]];if(!length(j))return(NULL)
    data.table(zip=zips$zip[i],SID=sw$SID[j],NAME=sw$NAME[j],Wind=sw$Wind[j],shared=as.IDate(sw$Date0[j]))
  }))
}
note("Assigning ZIP intersections at both radii")
original<-assign_pairs(swaths);large<-assign_pairs(expanded)
saved<-fread(protected[1],colClasses=c(zip="character"))
saved[,shared:=as.IDate(shared)]
keys<-c("zip","SID","NAME","Wind","shared")
stopifnot(fsetequal(original[,..keys],saved[,..keys]),
          !anyDuplicated(large[,.(zip,SID)]),
          nrow(original[!large,on=.(zip,SID)])==0L,
          fsetequal(unique(original[,.(SID,NAME,Wind,shared)]),unique(large[,.(SID,NAME,Wind,shared)])))
fwrite(large,file.path(out,"zip_storm_pairs_rmw_1p5.csv"))
fwrite(meta,file.path(out,"storm_inventory.csv"))
saveRDS(list(original=swaths,expanded=expanded),file.path(out,"swaths.rds"))
note("Preparing cached annual population for expanded ZIP universe")
files<-list.files("output/event_time",pattern="^acs_[0-9]+[.]rds$",full.names=TRUE)
annual<-rbindlist(lapply(files,function(f) {
  a<-as.data.table(readRDS(f))
  a[,.(zip=GEOID,year=as.integer(sub(".*acs_([0-9]+)[.]rds","\\1",f)),pop=estimate)]
}))[zip%in%large$zip]
population<-annual[,.(year=2005:2022,pop=approx(year,pop,xout=2005:2022,method="linear",rule=2)$y),by=zip]
stopifnot(!anyDuplicated(population[,.(zip,year)]),all(is.finite(population$pop)),
          all(unique(large$zip)%in%population$zip))
savedpop<-fread(protected[2],colClasses=c(zip="character"))
check<-merge(savedpop,population,by=c("zip","year"),suffixes=c("_old","_new"))
poperr<-max(abs(check$pop_old-check$pop_new))
stopifnot(nrow(check)==nrow(savedpop),poperr<1e-8)
fwrite(population,file.path(out,"population_annual.csv"))
audit<-rbindlist(lapply(c(34,55,64),function(t) {
  a<-original[Wind>=t & year(shared)>=2005];b<-large[Wind>=t & year(shared)>=2005]
  data.table(threshold=t,original_zips=uniqueN(a$zip),expanded_zips=uniqueN(b$zip),
    original_pairs=nrow(a),expanded_pairs=nrow(b),storms=uniqueN(a$SID),
    new_pairs=nrow(b[!a,on=.(zip,SID)]))
}))
fwrite(audit,file.path(out,"exposure_comparison.csv"))
writeLines(c("Original-radius ZIP assignments exactly reproduce the saved assignments.",
  "Every original ZIP–storm pair is included at 1.5 RMW.",
  "Storm dates, intensities, retained track segments and storm inventory are held fixed.",
  sprintf("Original-ZIP annual population maximum difference: %.3g",poperr),
  "Expanded swaths change geographic assignment only; no local timing analysis is performed."),
  file.path(out,"geometry_validation.txt"))
print(audit)
