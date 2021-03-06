library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggrepel)

nExamples <- 8

# load data
countriesFit <- readRDS("./generated_data/countriesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI) %>%
  as_tibble(.) %>% 
  mutate(Country=ifelse(Country %in% "Korea, South","South Korea",levels(Country)[Country]))
  

simFit <- readRDS("./generated_data/SIR_dynamics_fit.RDS") %>%
  dplyr::mutate(allee = factor(ifelse(allee,"with NPIs","without NPIs")))

fittedCountries <- as.character(countriesFit$Country)

countriesFit$allee = factor("Countries", levels = c("Countries", levels(simFit$allee)))
simFit$Id <- seq(10000, 10000 + nrow(simFit)-1, 1)
 
mergeCols <- c("Id", "allee", "slopeI", "slopeF", "slopeRatio")
#("Id","allee","slopeI","slopeF","slopeRatio")

# get countries with the least and most pronounced breakpoints
increasingRatio <- dplyr::filter(countriesFit, slopeRatio > 1)
slopeOrder <- base::order(increasingRatio$slopeRatio)
nCountries <- length(slopeOrder)
#ind <- round(c(seq(1, nCountries, nCountries/(nExamples-1)), nCountries))
ind <- c(seq(1, nExamples/2, 1), seq(nCountries-nExamples/2, nCountries, 1))
examplesIndices <- slopeOrder[ind]
examplesId <- as.character(increasingRatio$Country[examplesIndices])

countriesFit$label <- NA
countriesFit$Country <- as.character(countriesFit$Country)
countriesFit <- countriesFit %>%
  dplyr::mutate(., label = ifelse(Country %in% examplesId, as.character(Country), label)) %>%
  as_tibble(.)

sim_countriesFit <- merge(simFit[,mergeCols], countriesFit, by=mergeCols, all=T) %>%
  as_tibble(.)

sim_countriesFit$breakout <- sim_countriesFit$slopeRatio > 1
sim_countriesFit$breakout_allee = interaction(sim_countriesFit$allee, sim_countriesFit$breakout)

#### ----- PLOT FIGURA 2 ------

xlims=c(0,.39)
ylims= c(-.01,2)

plot_slopes <- sim_countriesFit %>%
  dplyr::mutate_at(vars(contains("slope"),"slopeRatio"), .funs = list("log10"=function(x) log10(1+x))) %>%
  ggscatterhist(x="slopeI_log10",y="slopeRatio_log10",
                #color="bajoUmbral",
                color="allee",
                fill="allee",
                shape="allee",
                size="allee",
#                alpha="allee",
                group="allee",
                alpha=.4,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(1+slope before breakout point)",
                ylab="log10(1+slope after/slope before)",
                palette = c("#b33018","#14b74b","#000080"),
                ylim=ylims,
                xlim=xlims)

plot_slopes$sp <- simFit %>%
  ggplot(., aes(group=allee, fill=allee, x=log10(1+slopeI), y=log10(1+slopeRatio))) +
  stat_density_2d(geom="polygon", aes(fill=allee, alpha=..level..), contour=T, bins=14) + 
  geom_point(shape=20, data=(countriesFit), alpha=1,
             aes(col=allee), size = 0.8) +
  scale_fill_manual(values = c("#FFFFFF","#b33018","#14b74b")) +
  theme_classic() +
  scale_color_manual(values=c("#000080","#b33018","#14b74b"))+
  theme(legend.position = "top") +  
  geom_abline(slope=0,intercept = log10(2),linetype="dotted") +
  geom_text_repel(data=sim_countriesFit, aes(label=label), size=3)+
  annotate("text", x=0.04, y=log10(2)-0.08, label="equal slopes", size=3) +
  xlim(xlims) +
  ylim(ylims) +
#  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
  labs(x=expression(paste(Log[10], "(", 1+mu[1], ")")),
       y=expression(paste(Log[10], "(", 1+frac(mu[2],mu[1]), ")"))) 


plot_slopes$sp$labels$colour=""
plot_slopes$sp$labels$fill=""
plot_slopes$sp$labels$shape=""
plot_slopes$sp$labels$size=""
plot_slopes$sp$labels$alpha=""
plot_slopes$yplot <- plot_slopes$yplot + ylim(ylims)
plot_slopes$xplot <- plot_slopes$xplot + ylim(xlims)
plot_slopes$sp <- plot_slopes$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes <- print(plot_slopes)

scatterCountriesPlot <- plot_slopes

#saveRDS(list(scatterSlopesCountries = plot_slopes), "./generated_data/countries_scatter_plot.RDS")
ggsave(plot_slopes, file="fig2_scatter_condados.pdf", width=5, height=4, units =  "in")
#ggsave(plot_slopes, file="fig2_scatter_condados.png", width=5, height=4)

