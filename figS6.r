# - FIGURA SUPPLEMENTARY MATERIALS S6 - Simulated allee dynamics: evidence for breaking point and \mu parameters distribution

library(ggpubr)
library(ggforce)
library(ggrepel)



#fitCoef - vien de Alle_dynamcis_simulation.R
figS6a = fitCoefs %>%#filter(weighted.evidence>.8) %>%  
  ggplot(aes(x=weighted.evidence,fill=allee.f))+
  geom_histogram() +
  facet_grid(~allee.f)+
  labs(x="Weighted evidence for breakpoint",y="# of simulations",fill="")+
  scale_fill_manual(values = c("#b33018", "#14b74b"))+
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size=10),
    text = element_text(size=12),
    legend.position = "top")  

figS6a




xlims=c(0,.45)
ylims= c(-.01,2.5)


plot_slopes_s6b=
  fitCoefs %>% 
  mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
  #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
  filter(cociente>1) %>% 
  ggscatterhist(x="slopeI_log10",y="cociente_log10",
                #color="bajoUmbral",
                color="allee.f",
                fill="allee.f",
                size="allee.f",
                #                alpha="allee.f",
                group="allee.f",
                alpha=.4,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(1+slope before breakout point)",
                ylab="log10(1+slope after/slope before)",
                palette = c("#14b74b","#b33018"),
                ylim=ylims,
                xlim=xlims
  )


plot_slopes_s6b$sp = 
  fitCoefs %>% filter(cociente>1) %>%
  ggplot(aes(group=allee.f,fill=allee.f,x=log10(1+slopeI),y=log10(1+cociente))) +
  stat_density_2d(geom="polygon",aes(fill=allee.f,alpha=..level..),contour=T,bins=20) + 
    scale_fill_manual(values = c("#14b74b","#b33018")) + theme_classic() +
  scale_color_manual(values=c("#14b74b","#b33018"))+
  theme(legend.position = "none") +  
  geom_abline(slope=0,intercept = log10(2)-.1,linetype="dotted") +
  xlim(xlims) + ylim(ylims) +
  #  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
  labs(x=expression(paste(Log10, "(", 1+mu[1], ")")),
       y=expression(paste(Log10, "(", 1+frac(mu[1],mu[2]), ")"))) 

plot_slopes_s6b$sp$labels$colour=""
plot_slopes_s6b$sp$labels$fill=""
plot_slopes_s6b$sp$labels$shape=""
plot_slopes_s6b$sp$labels$size=""
plot_slopes_s6b$sp$labels$alpha=""
plot_slopes_s6b$yplot=plot_slopes_s6b$yplot + ylim(ylims) + theme(legend.position = "none")
plot_slopes_s6b$xplot=plot_slopes_s6b$xplot + ylim(xlims)+ theme(legend.position = "none")
plot_slopes_s6b$sp = plot_slopes_s6b$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes_s6b=print(plot_slopes_s6b )

fig26 = ggarrange(figS6a,plot_slopes_s6b,legend = "top",labels = "auto")
fig26

ggsave(fig26,file="figs6.pdf",width=8,height=4)
ggsave(fig26,file="figs6.png",width=8,height=4)
