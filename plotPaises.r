bb_p=as_tibble(b)


colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")

load('cuadrantes_paises.RData')
bb_p$allee.f=factor("countries",levels = c(levels(cuadrantes_paises$allee.f),"countries"))
plot_bbp_slopes=
  bb_p %>% mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                color="allee.f",
                #color="#000080",
                alpha=.3,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
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
plot_bbp_slopes

plot_bbp_slopes = print(plot_bbp_slopes)


ggsave(plot_bbp_slopes,file="fig2_paises_slopes_2.pdf")
  

#histograma poblacion paises
bb_p %>% ggplot(aes(x=log(Population)))+geom_histogram()
bb_p %>% summarize(m=mean(log(Population),na.rm=T),s=sd(log(Population),na.rm=T))

