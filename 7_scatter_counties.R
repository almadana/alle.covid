library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggrepel)

# load data
countiesFit <- readRDS("./generated_data/countiesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI)

simFit <- readRDS("./generated_data/SIR_dynamics_fit.RDS") %>%
  dplyr::mutate(allee = factor(ifelse(allee,"with Allee effect","without Allee effect")))

fittedCounties <- as.character(countiesFit$Id)

stateAbrev <- read.csv('county_abbreviations.csv')

countiesCovid <- read.csv("us-counties.csv", colClasses = c("fips" = "character")) %>%
  as_tibble(.) %>%
  dplyr::filter(., county != "Unknown") %>%
  dplyr::rename(., Id = fips) %>%
  dplyr::group_by(Id) %>% 
  dplyr::mutate(t = row_number(), peak = which.max(diff(cases))[1]) %>% 
  dplyr::filter(t > 14, t <= peak) %>% 
  dplyr::mutate(t = t - 14) %>%
  merge(., stateAbrev, by = "state", all.x = T) %>%
  dplyr::mutate(., county.full = paste(county, code, sep=", ")) %>%
  dplyr::filter(., Id %in% fittedCounties) %>%
  tidyr::as_tibble(.)

county_short <- countiesCovid %>%
  group_by(Id) %>%
  summarize(county.full = county.full[1]) %>%
  tidyr::as_tibble(.)

countiesFit <- merge(countiesFit, county_short, by = "Id", all.x = T)

countiesFit$label <- NA
county4 <- NULL
countiesFit <- countiesFit %>%
  dplyr::mutate(label = ifelse(Id %in% county4, county.full, label)) %>%
  as_tibble(.)

countiesFit$allee = factor("U.S. counties", levels = c("U.S. counties", levels(simFit$allee)))
simFit$Id <- seq(10000, 10000 + nrow(simFit)-1, 1)
 
mergeCols <- c("Id", "allee", "slopeI", "slopeF", "slopeRatio")
#("Id","allee","slopeI","slopeF","slopeRatio")

sim_countiesFit <- merge(simFit[,mergeCols], countiesFit, by=mergeCols, all=T) %>%
  as_tibble(.)

sim_countiesFit$breakout <- sim_countiesFit$slopeRatio > 1
sim_countiesFit$breakout_allee = interaction(sim_countiesFit$allee, sim_countiesFit$breakout)

#### ----- PLOT FIGURA 2 ------

xlims=c(0,.35)
ylims= c(-.01,2)

plot_slopes <- sim_countiesFit %>%
  dplyr::mutate_at(vars(contains("slope"),"slopeRatio"), .funs = list("log10"=function(x) log10(1+x))) %>%
  #filtra a los que no tiene quiebre de pendientes positivo (ratio mayor a uno)
  #filter(slopeRatio>1) %>% 
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
  filter(slopeRatio>1) %>%
  ggplot(., aes(group=allee, fill=allee, x=log10(1+slopeI), y=log10(1+slopeRatio))) +
  stat_density_2d(geom="polygon", aes(fill=allee, alpha=..level..), contour=T, bins=14) + 
  geom_point(shape=20, data=(countiesFit %>% filter(slopeRatio>1)), alpha=1, aes(col=allee)) +
  scale_fill_manual(values = c("#FFFFFF","#b33018","#14b74b")) +
  theme_classic() +
  scale_color_manual(values=c("#000080","#b33018","#14b74b"))+
  theme(legend.position = "top") +  
  geom_abline(slope=0,intercept = log10(2)-.1,linetype="dotted") +
  geom_text_repel(data=sim_countiesFit, aes(label=label), size=3)+
  annotate("text", x=0.03, y=log10(2)+0.1, label="equal slopes", size=3) +
  xlim(xlims) +
  ylim(ylims) +
#  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
  labs(x=expression(paste(Log10, "(", 1+mu[1], ")")),
       y=expression(paste(Log10, "(", 1+frac(mu[2],mu[1]), ")"))) 


plot_slopes$sp$labels$colour=""
plot_slopes$sp$labels$fill=""
plot_slopes$sp$labels$shape=""
plot_slopes$sp$labels$size=""
plot_slopes$sp$labels$alpha=""
plot_slopes$yplot <- plot_slopes$yplot + ylim(ylims)
plot_slopes$xplot <- plot_slopes$xplot + ylim(xlims)
plot_slopes$sp <- plot_slopes$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes <- print(plot_slopes)

saveRDS(list(scatterSlopesCounties = plot_slopes), "./generated_data/counties_scatter_plot.RDS")
ggsave(plot_slopes, file="fig2_scatter_condados.pdf", width=5, height=4)
#ggsave(plot_slopes, file="fig2_scatter_condados.png", width=5, height=4)


