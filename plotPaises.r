bb_p=as_tibble(b)


colnames(cuadrantes)[c(8,10)]=c("initial.slope","slope.after.thresh")


plot_bbp_slopes=
  bb_p %>% mutate_at(c("initial.slope","slope.after.thresh"),.funs = list("log10"=function(x) log10(1+x))) %>% 
  ggscatterhist(x="initial.slope",y="slope.after.thresh",
                #color="bajoUmbral",
                color="#000080",
                alpha=.4,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="slope after threshold",
                #palette = "#000080",
                palette = c("#ff9e47","#14b74b","#000080"),
                
                ylim=c(0,10),
                xlim=c(0,2)
  )
plot_bbp_slopes$sp = plot_bbp_slopes$sp +
  stat_chull(data=cuadrantes,aes(fill=allee.f),alpha=0.5,geom="polygon") + 
  geom_abline(intercept = 0,slope = 1)
plot_bbp_slopes$sp$labels$colour=""
plot_bbp_slopes


#histograma poblacion paises
bb_p %>% ggplot(aes(x=log(Population)))+geom_histogram()
bb_p %>% summarize(m=mean(log(Population),na.rm=T),s=sd(log(Population),na.rm=T))

