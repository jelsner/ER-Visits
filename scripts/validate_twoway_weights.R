# Check the product-weight implementation against explicitly copied clusters.
suppressPackageStartupMessages({library(data.table);library(fixest)})
setDTthreads(1);setFixest_nthreads(1)
expr <- parse("scripts/regenerate_suppression_figures.R")
for(e in expr) if(is.call(e)&&identical(e[[1]],as.name("<-")) &&
  is.symbol(e[[2]])&&as.character(e[[2]])%in%c("aggregate_cells","fit_pair","gcomp")) eval(e)
set.seed(314159)
d <- CJ(zip=sprintf("%05d",1:20),season=2017:2019,day=1:42)
d[,year_month:=paste0(season,"-06")]
d[,dow:=(day-1L)%%7L+1L]
d[,event:=factor(rep(c("None","Threat","Impact","Cleanup"),
                      times=c(18,8,8,8))[day],
                  levels=c("None","Threat","Impact","Cleanup"))]
d[,pop:=1000+as.integer(zip)*40+(season-2017)*30]
d[,lambda:=7+as.integer(zip)*.15+fifelse(event=="Impact",-2,
              fifelse(event=="Cleanup",1,fifelse(event=="Threat",.5,0))) ]
d[,y:=rpois(.N,lambda)]
d[,`:=`(suppressed=y<=5,er_visits=fifelse(y>=6,as.numeric(y),NA_real_))]
zz <- sprintf("%05d",1:20); years<-as.character(2017:2019)
zip_draw <- sample(zz,length(zz),replace=TRUE)
season_draw <- c("2017","2017","2018")
wz <- table(factor(zip_draw,levels=zz));wy<-table(factor(season_draw,levels=years))
product <- fit_pair(aggregate_cells(d),wz,wy)
explicit <- rbindlist(lapply(seq_along(zip_draw),function(i)
  rbindlist(lapply(seq_along(season_draw),function(j) {
    q<-copy(d[zip==zip_draw[i] & season==as.integer(season_draw[j])])
    q[,`:=`(zip=sprintf("copy_%02d",i),year_month=paste0("season_copy_",j))];q
  }))))
copied <- fit_pair(aggregate_cells(explicit))
levels <- c("Threat","Impact","Cleanup")
a <- gcomp(product,levels);b<-gcomp(copied,levels)
coef_err <- max(abs(coef(product$a)-coef(copied$a)),
                abs(coef(product$b)-coef(copied$b)))
contrast_err <- max(abs((a$C+3.8*a$A)-(b$C+3.8*b$A)))
stopifnot(coef_err<1e-6,contrast_err<1e-6,
          sum(product$target$n*product$target$mult)==nrow(d),
          nrow(explicit)==nrow(d))
result <- data.table(max_coefficient_error=coef_err,max_contrast_error=contrast_err,
                     weighted_rows=sum(product$target$n*product$target$mult),
                     explicit_rows=nrow(explicit))
fwrite(result,"output/twoway_bootstrap/product_weight_validation.csv")
print(result)
