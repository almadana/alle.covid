bb_p=as_tibble(b)



load('cuadrantes_paises.RData')
bb_p$allee.f=factor("Countries & regions",levels = c(levels(cuadrantes_paises$allee.f),"Countries & regions"))

xlims=c(-0.25,1)
ylims= c(-0,1)


bb_p$slope.after.thresh = bb_p$initial.slope+ bb_p$slope.after.thresh
bb_p$cociente = bb_p$slope.after.thresh / bb_p$initial.slope

#cuadrantes_paises$slope.after.thresh =cuadrantes_paises$initial.slope +  cuadrantes_paises$slope.after.thresh
bb_p$Id


country_names = country_data %>% group_by(Id) %>% summarize(country=country[1])
bb_p=merge(bb_p,country_names,by="Id")

cuadrantes_paises$Id=seq(10000,10000+nrow(cuadrantes_paises)-1,1)

sim_bb_p = merge(cuadrantes_paises[,c("Id","n","allee.f","initial.slope","slope.after.thresh","cociente")],bb_p,by=c("Id","allee.f","initial.slope","slope.after.thresh","cociente"),all=T)
View(sim_bb_p)
sim_bb_p$allee.f=relevel(sim_bb_p$allee.f,"Countries & regions")
sim_bb_p$breakout = sim_bb_p$cociente>1
sim_bb_p$breakout_allee.f = interaction(sim_bb_p$allee.f,sim_bb_p$breakout)

sim_bb_p$label=NA
sim_bb_p=sim_bb_p %>% mutate(label=ifelse(Id %in% country4,country,label))


xlims=c(0,.5)
ylims= c(-.01,3)



plot_slopes_p=
  sim_bb_p %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
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
                palette = c("#4020ab","#b33018","#14b74b"),
                ylim=ylims,
                xlim=xlims
  )


plot_slopes_p$sp = 
  cuadrantes_paises %>% filter(cociente>1) %>%
  ggplot(aes(group=allee.f,fill=allee.f,x=log10(1+initial.slope),y=log10(1+cociente))) +
  stat_density_2d(geom="polygon",aes(fill=allee.f,alpha=..level..),contour=T,bins=14) + 
  geom_point(shape=20,data=(bb_p%>% filter(cociente>1)),alpha=1,aes(col=allee.f))+
  scale_fill_manual(values = c("#FFFFFF","#b33018","#14b74b")) + theme_classic() +
  scale_color_manual(values=c("#4020ab","#b33018","#14b74b"))+
  theme(legend.position = "top") +  
  geom_abline(slope=0,intercept = log10(2),linetype="dotted") +
  annotate("text",x=0.03,y=log10(2)+0.1,label="equal slopes",size=3) +
  xlim(xlims) + ylim(ylims) +
  geom_text_repel(data=sim_bb_p,aes(label=label),size=3)+
  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")


plot_slopes_p$sp$labels$colour=""
plot_slopes_p$sp$labels$fill=""
plot_slopes_p$sp$labels$shape=""
plot_slopes_p$sp$labels$size=""
plot_slopes_p$sp$labels$alpha=""
plot_slopes_p$yplot=plot_slopes_p$yplot + ylim(ylims)
plot_slopes_p$xplot=plot_slopes_p$xplot + ylim(xlims)
plot_slopes_p$sp = plot_slopes_p$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes_p=print(plot_slopes_p)

ggsave(plot_slopes_p,file="fig2_scatter_paises.pdf",width=5,height=4)
ggsave(plot_slopes_p,file="fig2_scatter_paisess.png",width=5,height=4)


# 
# plot_bbp_slopes=
#   bb_p %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
#   filter(cociente_log10<10) %>% 
#   ggscatterhist(x="initial.slope_log10",y="cociente_log10",
#                 color="allee.f",
#                 #color="#000080",
#                 alpha=.3,size=3,
#                 margin.plot="boxplot",
#                 ggtheme=theme_bw(),
#                 xlab="initial slope",
#                 ylab="slope after threshold",
# 
#                 palette = "#000080",
#                 #palette = c("#b33018","#14b74b","#000080"),
#                 
#                 #                ylim=c(0,15),
#                 #               xlim=c(0,2)
#   )
# plot_bbp_slopes$sp = plot_bbp_slopes$sp  +  stat_chull(data=cuadrantes_paises,aes(fill=allee.f),alpha=.3,geom="polygon") + 
#   geom_abline(intercept = 0,slope = 1) +scale_fill_manual(values = c("#b33018","#14b74b","#000080"))
# plot_bbp_slopes$sp$labels$colour=""
# plot_bbp_slopes$sp$labels$fill="simulations"
# plot_bbp_slopes$yplot = plot_bbp_slopes$yplot + ylim(c(0,12))
# =======
#                 #palette = "#000080",
#                 palette = c("#0040dd"),
#                 ylim=ylims,
#                 xlim=xlims
#  #               ylim=c(0,15),
# #                xlim=c(0,1.5)
#   )
# plot_bbp_slopes$yplot = plot_bbp_slopes$yplot + ylim(ylims)
# plot_bbp_slopes$xplot = plot_bbp_slopes$xplot + ylim(xlims)
# 
# plot_bbp_slopes$sp = plot_bbp_slopes$sp +
#   geom_mark_hull(concavity = 5,radius=.035,aes(fill=allee.f)) + 
#   geom_abline(intercept = 0,slope = 1,linetype="dotted")
# plot_bbp_slopes$sp$labels$colour=""
# plot_bbp_slopes$sp$labels$fill=""
# 
# >>>>>>> alvaro
# plot_bbp_slopes
# plot_bbp_slopes = print(plot_bbp_slopes)
# 
# <<<<<<< HEAD
# plot_bbp_slopes = print(plot_bbp_slopes)
# 
# 
# ggsave(plot_bbp_slopes,file="fig2_paises_slopes_2.pdf")
#   
# =======
# ggsave(plot_bbp_slopes,file="fig2_paises.pdf",width = 5,height=4)
# >>>>>>> alvaro
# 
# #histograma poblacion paises
# bb_p %>% ggplot(aes(x=log(Population)))+geom_histogram()
# bb_p %>% summarize(m=mean(log(Population),na.rm=T),s=sd(log(Population),na.rm=T))
# 
