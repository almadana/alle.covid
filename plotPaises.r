bb_p=as_tibble(b)


colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")

<<<<<<< HEAD
load('cuadrantes_paises.RData')
bb_p$allee.f=factor("countries",levels = c(levels(cuadrantes_paises$allee.f),"countries"))
=======
bb_p$allee.f=factor("Countries",levels = c(levels(cuadrantes$allee.f),"Countries"))

xlims=c(-0.25,2)
ylims= c(-1,12)

>>>>>>> alvaro
plot_bbp_slopes=
  bb_p %>% mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                color="allee.f",
<<<<<<< HEAD
                #color="#000080",
                alpha=.3,size=3,
=======
                #color="#202080",
                alpha=.4,size=3,
>>>>>>> alvaro
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
<<<<<<< HEAD
                palette = "#000080",
                #palette = c("#b33018","#14b74b","#000080"),
                
                #                ylim=c(0,15),
                #               xlim=c(0,2)
  )
plot_bbp_slopes$sp = plot_bbp_slopes$sp  +  stat_chull(data=cuadrantes_paises,aes(fill=allee.f),alpha=.3,geom="polygon") + 
  geom_abline(intercept = 0,slope = 1) +scale_fill_manual(values = c("#b33018","#14b74b","#000080"))
plot_bbp_slopes$sp$labels$colour=""
plot_bbp_slopes$sp$labels$fill="simulations"
plot_bbp_slopes$yplot = plot_bbp_slopes$yplot + ylim(c(0,12))
=======
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

>>>>>>> alvaro
plot_bbp_slopes
plot_bbp_slopes = print(plot_bbp_slopes)

<<<<<<< HEAD
plot_bbp_slopes = print(plot_bbp_slopes)


ggsave(plot_bbp_slopes,file="fig2_paises_slopes_2.pdf")
  
=======
ggsave(plot_bbp_slopes,file="fig2_paises.pdf",width = 5,height=4)
>>>>>>> alvaro

#histograma poblacion paises
bb_p %>% ggplot(aes(x=log(Population)))+geom_histogram()
bb_p %>% summarize(m=mean(log(Population),na.rm=T),s=sd(log(Population),na.rm=T))

