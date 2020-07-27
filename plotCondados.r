library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggrepel)

bb=as_tibble(b)
#   
# bb$bajoUmbral = factor(bb$I.active<10)
# 
# bb %>% 
#   ggplot(aes(x=slope.after.thresh,fill=bajoUmbral))+geom_histogram(bins = 20,position = "identity",alpha=.5)
# 
# bb %>% 
#   ggplot(aes(x=initial.slope,fill=bajoUmbral))+geom_histogram(bins = 20,position = "identity",alpha=.5)
# 
# bb %>% group_by(bajoUmbral) %>% summarize(m=mean(slope.after.thresh),mx=max(slope.after.thresh))

bb$slope.after.thresh = bb$slope.after.thresh +bb$initial.slope
bb$cociente = bb$slope.after.thresh / bb$initial.slope

county_short= county_data %>% group_by(Id) %>% summarize(county.full=county.full[1])

bb = merge(bb,county_short,by="Id",all.x = T)

bb$label=NA
bb=bb %>% mutate(label=ifelse(Id %in% county4,county.full,label))


# #datos simulados
# load('simulados_allee.RData')


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


bb$allee.f=factor("U.S. counties",levels = c("U.S. counties",levels(cuadrantes$allee.f)))
# 
# xlims=c(0,.5)
# ylims= c(0,2)
# plot_bb_slopes=
#   bb %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>% 
#     ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 color="allee.f",
#                 #color="#000080",
#                 alpha=.4,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="initial slope",
#                 ylab="slope after threshold",
#                 #palette = "#000080",
#                 palette = c("#000080"),
#                 
#                 ylim=ylims,
#                 xlim=xlims
#   )
# # plot_bb_slopes$sp = plot_bb_slopes$sp +
# #   stat_chull(data=cuadrantes_paises,aes(fill=allee.f),alpha=0.5,geom="polygon") + 
# #   geom_abline(intercept = 0,slope = 1)
# plot_bb_slopes$sp = 
#   plot_bb_slopes$sp +
#     geom_mark_hull(concavity = 2,aes(fill=allee.f)) + 
# #  stat_chull(aes(fill=allee.f),alpha=0.5,geom="ellipse") + 
#     geom_abline(intercept = log10(2),slope = 0,linetype="dotted") 
# 
# plot_bb_slopes$sp$labels$colour=""
# plot_bb_slopes$sp$labels$fill=""
# plot_bb_slopes$yplot = plot_bb_slopes$yplot + ylim(ylims)
# plot_bb_slopes$xplot = plot_bb_slopes$xplot + ylim(xlims) # OJO!! ES UN BARPLOT ROTADO 90° - PEÑAROL INTELIGENCIA
# 
# plot_bb_slopes
# 
# plot_bb_slopes =print(plot_bb_slopes)
# 
# ggsave(plot_bb_slopes,file="fig2_counties.pdf",width = 5,height = 4)
# 
# colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")
# colnames(bb)
# cuadrantes$Id=seq(10000,10000+nrow(cuadrantes)-1,1)
# bb$allee.f="U.S. counties"
# 

sim_bb = merge(cuadrantes[,c("Id","n","allee.f","initial.slope","slope.after.thresh","cociente")],bb,by=c("Id","allee.f","initial.slope","slope.after.thresh","cociente"),all=T)
View(sim_bb)
sim_bb$allee.f=relevel(sim_bb$allee.f,"U.S. counties")
sim_bb$breakout = sim_bb$cociente>1
sim_bb$breakout_allee.f = interaction(sim_bb$allee.f,sim_bb$breakout)

#### ----- PLOT FIGURA 2 ------

xlims=c(0,.35)
ylims= c(-.01,2)



plot_slopes=
  sim_bb %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
  #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
  filter(cociente>1) %>% 
  ggscatterhist(x="initial.slope_log10",y="cociente_log10",
                #color="bajoUmbral",
                color="allee.f",
                fill="allee.f",
                shape="allee.f",
                size="allee.f",
#                alpha="allee.f",
                group="allee.f",
                alpha=.4,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(1+slope before breakout point)",
                ylab="log10(1+slope after/slope before)",
                palette = c("#000080","#b33018","#14b74b"),
                ylim=ylims,
                xlim=xlims
  )

# plot_slopes$sp = plot_slopes$sp + scale_shape_manual(values=c(4,20,20)) +
#   scale_size_manual(values=c(2,2,2))+#scale_alpha_manual(values=c(.5,.3,.3)) +
#      geom_abline(slope=0,intercept = log10(2),linetype="dotted") +
#   annotate("text",x=0.03,y=log10(2)+0.1,label="equal slopes",size=3)+theme_classic() +
#   theme(legend.position = "top")
#   
# 
# plot_slopes_solo_condados=
#   sim_bb %>% filter(allee.f=="U.S. counties") %>% 
#   mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
#   filter(cociente>1) %>% 
#   ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 #color="bajoUmbral",
#                 color="allee.f",
#                 alpha=.4,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="log10(1+initial slope)",
#                 ylab="log10(1+slope after threshold/initial slope)",
#                 palette = c("#000080"),
#                 ylim=ylims,
#                 xlim=xlims
#   )
# 

# plot con división en el cociente
# 
# plot_slopes_breakout=
#   sim_bb %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
#   # filter(cociente>1) %>% 
#   ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 #color="bajoUmbral",
#                 color="breakout_allee.f",
#                 fill="allee.f",
#                 shape="allee.f",
#                 size="allee.f",
#                 alpha="allee.f",
#                 alpha=.4,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="log10(1+initial slope)",
#                 ylab="log10(1+slope after threshold/initial slope)",
#                 palette = c("#000080","#b33018","#14b74b","#000080","#b33018","#14b74b"),
#                 ylim=ylims,
#                 xlim=xlims
#   )
# 
# #plot_slopes_breakout$yplot = plot_slopes_breakout$yplot + ylim(ylims)
# 
# plot_slopes$yplot = plot_slopes_breakout$yplot

# plot_slopes_con_allee=
#   sim_bb %>% filter(allee.f=="with Allee effect") %>% 
#   mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
#   filter(cociente>1) %>% 
#   ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 #color="bajoUmbral",
#                 color="allee.f",
#                 alpha=.4,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="log10(1+initial slope)",
#                 ylab="log10(1+slope after threshold/initial slope)",
#                 palette = c("#b33018"),
#                 ylim=ylims,
#                 xlim=xlims
#   )
# 
# plot_slopes_con_allee$sp
# 
# plot_slopes_sin_allee=
#   sim_bb %>% filter(allee.f=="without Allee effect") %>% 
#   mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
#   filter(cociente>1) %>% 
#   ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 #color="bajoUmbral",
#                 color="allee.f",
#                 alpha=.4,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="log10(1+initial slope)",
#                 ylab="log10(1+slope after threshold/initial slope)",
#                 palette = c("#14b74b"),
#                 ylim=ylims,
#                 xlim=xlims
#   )
# 
# plot_slopes_sin_allee$sp
# 
# ggarrange(print(plot_slopes),print(plot_slopes_con_allee$sp+legend("null")),print(plot_slopes_sin_allee$sp),
#           nrow = 1,common.legend = T,legend="bottom")
# 
# plot_slopes$sp <-  plot_slopes_solo_condados$sp+
#   geom_abline(slope=0,intercept = log10(2),linetype="dotted") #+  
# # geom_mark_hull(concavity = 5,expand = -.00000001,aes(fill=allee.f)) 



#plot_slopes$sp + stat_density_2d(geom="polygon",aes(fill=allee.f,alpha=..level..),contour=T,bins=14)
  

plot_slopes$sp = 
  cuadrantes %>% filter(cociente>1) %>%
  ggplot(aes(group=allee.f,fill=allee.f,x=log10(1+initial.slope),y=log10(1+cociente))) +
  stat_density_2d(geom="polygon",aes(fill=allee.f,alpha=..level..),contour=T,bins=14) + 
  geom_point(shape=20,data=(bb%>% filter(cociente>1)),alpha=1,aes(col=allee.f))+
  scale_fill_manual(values = c("#FFFFFF","#b33018","#14b74b")) + theme_classic() +
  scale_color_manual(values=c("#000080","#b33018","#14b74b"))+
  theme(legend.position = "top") +  
  geom_abline(slope=0,intercept = log10(2)-.1,linetype="dotted") +
  geom_text_repel(data=sim_bb,aes(label=label),size=3)+
  annotate("text",x=0.03,y=log10(2)+0.1,label="equal slopes",size=3) +
  xlim(xlims) + ylim(ylims) +
#  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
  labs(x=expression(paste(Log10, "(", 1+mu[1], ")")),
       y=expression(paste(Log10, "(", 1+frac(mu[1],mu[2]), ")"))) 


plot_slopes$sp$labels$colour=""
plot_slopes$sp$labels$fill=""
plot_slopes$sp$labels$shape=""
plot_slopes$sp$labels$size=""
plot_slopes$sp$labels$alpha=""
plot_slopes$yplot=plot_slopes$yplot + ylim(ylims)
plot_slopes$xplot=plot_slopes$xplot + ylim(xlims)
plot_slopes$sp = plot_slopes$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes=print(plot_slopes)

ggsave(plot_slopes,file="fig2_scatter_condados.pdf",width=5,height=4)
ggsave(plot_slopes,file="fig2_scatter_condados.png",width=5,height=4)

# 
# plot_ch = cuadrantes %>% 
#   filter(cociente<10,cociente>log10(2)) %>% 
#   ggplot(aes(x=log10(1+initial.slope),y=log10(1+cociente),fill=allee.f,col=allee.f)) + geom_mark_hull(concavity = 5) +
#   #xlim(xlims)+ ylim(ylims)+
#   geom_point(size=2)
# plot_ch
#  