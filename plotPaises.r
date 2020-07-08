bb_p=as_tibble(b)


colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")

bb_p$allee.f=factor("Countries",levels = c(levels(cuadrantes$allee.f),"Countries"))

xlims=c(-0.25,2)
ylims= c(-1,12)

plot_bbp_slopes=
  bb_p %>% mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                color="allee.f",
                #color="#202080",
                alpha=.4,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
                #palette = "#000080",
                palette = c("#0040dd"),
                ylim=ylims,
                xlim=xlims
 #               ylim=c(0,15),
#                xlim=c(0,1.5)
  )
plot_bbp_slopes$yplot = plot_bbp_slopes$yplot + ylim(ylims)
plot_bbp_slopes$xplot = plot_bbp_slopes$xplot + ylim(xlims)

plot_bbp_slopes$sp = plot_bbp_slopes$sp +
  geom_mark_hull(concavity = 5,radius=.035,aes(fill=allee.f)) + 
  geom_abline(intercept = 0,slope = 1,linetype="dotted")
plot_bbp_slopes$sp$labels$colour=""
plot_bbp_slopes$sp$labels$fill=""

plot_bbp_slopes
plot_bbp_slopes = print(plot_bbp_slopes)

ggsave(plot_bbp_slopes,file="fig2_paises.pdf",width = 5,height=4)

#histograma poblacion paises
bb_p %>% ggplot(aes(x=log(Population)))+geom_histogram()
bb_p %>% summarize(m=mean(log(Population),na.rm=T),s=sd(log(Population),na.rm=T))

