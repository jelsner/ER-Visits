# Re-render the saved diagnostic without refitting the models.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
res <- fread("output/event_time/daily_estimates.csv")
summary <- fread("output/event_time/timing_summary.csv")
audit <- fread("output/event_time/sample_audit.csv")
log_note <- function(...) cat(...,"\n")
res[,threshold_label:=factor(paste0("\u2265",threshold," kt"),
                              levels=paste0("\u2265",c(34,55,64)," kt"))]
res[,alignment_label:=factor(alignment,levels=c("shared","local"),
             labels=c("Shared storm date (posted code)","ZIP-specific closest approach"))]
p <- ggplot(res,aes(day,difference,color=alignment_label,shape=alignment_label))+
  annotate("rect",xmin=-.5,xmax=.5,ymin=-Inf,ymax=Inf,fill="#f2e7d4",alpha=.7)+
  geom_hline(yintercept=0,color="#8b949e",linewidth=.4)+
  geom_vline(xintercept=0,color="#999999",linewidth=.4,linetype="dotted")+
  geom_line(linewidth=.8)+geom_point(size=2.1)+
  facet_wrap(~threshold_label,nrow=1)+
  scale_color_manual(values=c("#1864ab","#c05621"))+
  scale_shape_manual(values=c(16,17))+
  scale_x_continuous(breaks=seq(-8,8,2))+
  labs(title="Does the ED utilization trough align with impact day?",
       subtitle="Daily two-part model estimates within the same RMW storm-core swaths, 2005\u20132022",
       x="Days relative to the specified impact date",
       y="Difference in expected visits per ZIP-day",
       color=NULL,shape=NULL,
       caption=paste("Reference: ZIP-days outside all qualifying \u00b18-day windows. Population-weighted; common prediction support within each panel.",
                     "Point estimates only; no confidence intervals. Existing suppression and missing-row assumptions retained. Dates use UTC.",sep="\n"))+
  theme_minimal(base_size=12)+
  theme(legend.position="bottom",panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=17),
        strip.text=element_text(face="bold",size=13),
        plot.caption=element_text(hjust=0,size=9),plot.margin=margin(15,18,15,15))
ggsave("figs/daily_event_time_alignment.png",p,width=12.5,height=5.9,dpi=300,
       device=ragg::agg_png,bg="white")
si <- copy(res)
si[,alignment_label := factor(alignment,levels=c("shared","local"),
             labels=c("Shared storm date","ZIP-specific closest approach"))]
p_si <- p %+% si
p_si <- p_si + labs(title="Figure S1. Daily ED utilization by impact-date alignment")
ggsave("figs/FigS1_daily_event_time_alignment.png",p_si,
       width=12.5,height=5.9,dpi=300,device=ragg::agg_png,bg="white")
print(summary); print(audit)
capture.output(sessionInfo(),file="output/event_time/session_info.txt")
log_note("Saved daily event-time plot and estimates.")
