library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
bb=as_tibble(b)

bb$bajoUmbral = factor(bb$I.active<10)

bb %>% 
  ggplot(aes(x=slope.after.thresh,fill=bajoUmbral))+geom_histogram(bins = 20,position = "identity",alpha=.5)

bb %>% 
  ggplot(aes(x=initial.slope,fill=bajoUmbral))+geom_histogram(bins = 20,position = "identity",alpha=.5)

bb %>% group_by(bajoUmbral) %>% summarize(m=mean(slope.after.thresh),mx=max(slope.after.thresh))


#datos simulados
#load('simulados_allee.RData')


# plot_slopes = 
#   cuadrantes %>% mutate_at(c("coef.2","coef.4"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   ggscatterhist(x="coef.2",y="coef.4",
#                 color="allee.f",alpha=.3,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="initial slope",
#                 ylab="slope after threshold",
#                 palette = c("#ff9e47","#14b74b"),
#                 ylim=c(0,15),
#                 xlim=c(0,2)
#   )
# plot_slopes$sp <-plot_slopes$sp+geom_abline(slope=1,intercept = 0)
# plot_slopes$sp$labels$colour=""
# plot_slopes
# 
# 

bb$allee.f=factor("U.S. counties",levels = c(levels(cuadrantes$allee.f),"U.S. counties"))

xlims=c(-0.25,2.5)
ylims= c(-.5,8)
plot_bb_slopes=
bb %>% mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                color="allee.f",
                #color="#000080",
                alpha=.4,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
                #palette = "#000080",
                palette = c("#000080"),
                
                ylim=ylims,
                xlim=xlims
  )
# plot_bb_slopes$sp = plot_bb_slopes$sp +
#   stat_chull(data=cuadrantes_paises,aes(fill=allee.f),alpha=0.5,geom="polygon") + 
#   geom_abline(intercept = 0,slope = 1)
plot_bb_slopes$sp = 
  plot_bb_slopes$sp +
    geom_mark_hull(concavity = 5,radius=.035,aes(fill=allee.f)) + 
#  stat_chull(aes(fill=allee.f),alpha=0.5,geom="ellipse") + 
    geom_abline(intercept = 0,slope = 1,linetype="dotted") 

plot_bb_slopes$sp$labels$colour=""
plot_bb_slopes$sp$labels$fill=""
plot_bb_slopes$yplot = plot_bb_slopes$yplot + ylim(ylims)
plot_bb_slopes$xplot = plot_bb_slopes$xplot + ylim(xlims) # OJO!! ES UN BARPLOT ROTADO 90° - PEÑAROL INTELIGENCIA

plot_bb_slopes

plot_bb_slopes =print(plot_bb_slopes)

ggsave(plot_bb_slopes,file="fig2_counties.pdf",width = 5,height = 4)

colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")
colnames(bb)
cuadrantes$Id=seq(10000,10000+nrow(cuadrantes)-1,1)
bb$allee.f="counties"
sim_bb = merge(cuadrantes[,c("Id","n","allee.f","initial.slope","slope.after.thresh")],bb,by=c("Id","allee.f","initial.slope","slope.after.thresh"),all=T)
View(sim_bb)

plot_slopes=
  sim_bb %>% #mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                #color="bajoUmbral",
                color="allee.f",
                alpha=.4,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
                palette = c("#ff9e47","#14b74b","#000080"),
                ylim=c(0,15),
                xlim=c(0,2)
  )
plot_slopes$sp <-plot_slopes$sp+geom_abline(slope=1,intercept = 0) + 
  stat_chull(data=cuadrantes,aes(fill=allee.f),alpha=0.1,geom="polygon")
plot_slopes$sp$labels$colour=""

plot_slopes


plot_ch = cuadrantes %>% 
  ggplot(aes(x=initial.slope,y=slope.after.thresh,fill=allee.f)) + geom_mark_hull(concavity = 5)
plot_ch
