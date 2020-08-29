library(ggpubr)
library(ggforce)
library(ggrepel)
library(dplyr)

# load data
simFit <- readRDS("./generated_data/SIR_dynamics_fit.RDS") %>%
  dplyr::mutate(allee = factor(ifelse(allee,"with NPI","without NPI")))

countriesFit <- readRDS("./generated_data/countriesFit.RDS") %>%
  dplyr::as_tibble(.) 

countiesFit <- readRDS("./generated_data/countiesFit.RDS") %>%
  dplyr::as_tibble(.) 

dataFits <- list(counties = countiesFit, countries = countriesFit)


####################
# S4 make plots of distribution of slopes and correlation with population
####################

for (dataType in names(dataFits)) {
  fitData <- dataFits[[dataType]]
  fileName <- paste("./plots/S4_", dataType, "_fit_properties.png", sep = "")

  slopeI_hist <- ggplot(fitData, aes(slopeI)) +
    geom_histogram(color = "gray", fill="#4020ab") +
    xlab("Initial slope") +
    ylab("Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid"))
  #mtext("a",side=3,at=-2.2,cex=1.6)

  slopeF_hist <- ggplot(fitData, aes(slopeF)) +
    geom_histogram(color = "gray", fill="#4020ab") +
    xlab("Slope after threshold") +
    ylab("Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid"))
  #mtext("c",side=3,at=-2.2,cex=1.6)

  breakPointInfected_hist <- ggplot(fitData, aes(log10(I.active))) +
    geom_histogram(color = "gray", fill="#4020ab") +
    xlab("Log10(Active at threshold)") +
    ylab("Count") +
    #scale_x_continuous(trans = "log10") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid"))
  #mtext("b",side=3,at=-2.2,cex=1.6)

  slopeI_Pop_Corr <- cor.test(fitData$slopeI, log10(fitData$Population),
                             method = "pearson")
  slopeI_Pop_plot <- ggplot(fitData, aes(x = log10(Population), y = slopeI)) +
    geom_point(color = "#4020ab") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid")) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    ylab("Initial slope") +
    xlab("Log10(Population)")


  slopeF_Pop_Corr <- cor.test(fitData$slopeF, log10(fitData$Population),
                             method = "pearson")
  slopeF_Pop_plot <- ggplot(fitData, aes(x = log10(Population), y = slopeF)) +
    geom_point(color = "#4020ab") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid")) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    ylab("Slope after threshold") +
    xlab("Log10(Population)")


  breakPoint_Pop_Corr <- cor.test(log10(fitData$I.active), log10(fitData$Population),
                             method = "pearson")
  breakPoint_Pop_plot <- ggplot(fitData, aes(x = log10(Population),
                                                  y = log10(I.active))) +
    geom_point(color = "#4020ab") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid")) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    xlab("Log10(Active at threshold)") +
    xlab("Log10(Population)")

    
  figureSupp4 <- ggarrange(slopeI_hist, slopeF_hist, breakPointInfected_hist,
                   slopeI_Pop_plot, slopeF_Pop_plot, breakPoint_Pop_plot,
                   widths = c(4,3), labels = "auto")

  ggsave(figureSupp4, file = fileName, width = 11, height = 7, units = "in")
}


################
# S6 Histogram of weight of evidence of breaks for Allee and non-allee simulations
###############

# Plot histogram of the weight of evidence for breakpoint on simulations
# with and without allee
figS6a <- simFit %>%
  dplyr::filter(slopeRatio > 1) %>%
  ggplot(., aes(x = weighted.evidence, fill = allee))+
  geom_histogram() +
  facet_grid(~allee)+
  labs(x="Weighted evidence for breakpoint",y="# of simulations",fill="")+
  scale_fill_manual(values = c("#b33018", "#14b74b"))+
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size=10),
    text = element_text(size=12),
    legend.position = "top")  


#### plot density plot of change of slope as a function of initial slope.
#### a plot of log10(1 + slopeRatio) vs log10(1 + slopeI) is done

xlims=c(0,.45)
ylims= c(-.01,2)

plot_slopes_s6b <- simFit %>% 
  dplyr::mutate_at(., vars(contains("slope"), "slopeRatio"),
                   .funs = list("log10"=function(x) log10(1+x))) %>%
  # keep only those with positive slope ratio
  dplyr::filter(., !(is.nan(slopeRatio) | is.na(slopeRatio))) %>%
#  dplyr::filter(slopeRatio > 1) %>% 
  ggscatterhist(x="slopeI_log10", y="slopeRatio_log10",
                #color="bajoUmbral",
                color="allee",
                fill="allee",
                size="allee",
                #                alpha="allee",
                group="allee",
                alpha=.4,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(1+slope before breakout point)",
                ylab="log10(1+slope after/slope before)",
                palette = c("#b33018","#14b74b"),
                ylim=ylims,
                xlim=xlims)


plot_slopes_s6b$sp <- simFit %>%
  #dplyr::filter(slopeRatio>1) %>%
  ggplot(., aes(group = allee, fill = allee, x = log10(1+slopeI), y = log10(1+slopeRatio))) +
  stat_density_2d(geom="polygon", aes(fill=allee, alpha=..level..), contour=T, bins=14) + 
  scale_fill_manual(values = c("#b33018","#14b74b")) +
  theme_classic() +
  scale_color_manual(values=c("#b33018","#14b74b"))+
  theme(legend.position = "none") +  
  geom_abline(slope=0, intercept = log10(2), linetype="dotted") +
  xlim(xlims) +
  ylim(ylims) +
  #  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
  labs(x=expression(paste(Log10, "(", 1+mu[1], ")")),
       y=expression(paste(Log10, "(", 1+frac(mu[2],mu[1]), ")"))) 

plot_slopes_s6b$sp$labels$colour=""
plot_slopes_s6b$sp$labels$fill=""
plot_slopes_s6b$sp$labels$shape=""
plot_slopes_s6b$sp$labels$size=""
plot_slopes_s6b$sp$labels$alpha=""
plot_slopes_s6b$yplot=plot_slopes_s6b$yplot + ylim(ylims) + theme(legend.position = "none")
plot_slopes_s6b$xplot=plot_slopes_s6b$xplot + ylim(xlims)+ theme(legend.position = "none")
plot_slopes_s6b$sp = plot_slopes_s6b$sp + guides(alpha="none",shape="none",colour="none") 

plot_slopes_s6b <- print(plot_slopes_s6b)


#### put histograms and density plots together
figS6 <- ggarrange(figS6a, plot_slopes_s6b,
                   legend = "top",labels = "auto")

ggsave(figS6, file = "./plots/S6_weight_evidence_sims.png", width = 8, height = 4, units = "in")


################
# S? Histogram of weight of evidence of breaks for counties and countries
###############

weight.evid.vec <- c(countiesFit$weighted.evidence, countriesFit$weighted.evidence)
geoVec <- c(rep("U.S. counties", nrow(countiesFit)),
            rep("Countries & regions", nrow(countriesFit)))
dataEvidPlot <- list()

for (dataType in names(dataFits)) {
  fitData <- dataFits[[dataType]]
  fileName <- paste("./plots/SX_", dataType, "_weight_of_evidence.png", sep = "")

  if (dataType == "counties") {
    geoType <- "U.S. counties"  
  } else if (dataType == "countries") {
    geoType <- "Countries"  
  }
  fitData$geo <- geoType

  # Plot histogram of the weight of evidence for breakpoint on data
  dataEvidPlot[[dataType]] <- fitData %>%
    ggplot(., aes(x = weighted.evidence)) +
    geom_histogram() +
    labs(x="Weighted evidence for breakpoint",y="# of simulations",fill="")+
    scale_fill_manual(values = c("#b33018", "#14b74b"))+
    theme_classic() +
    theme(strip.background = element_blank(), strip.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5,size=10),
      text = element_text(size=12),
      legend.position = "top")  
}

figWeightEvid <- ggarrange(dataEvidPlot$counties, dataEvidPlot$countries,
                           labels = "auto", nrow = 1)

ggsave(figWeightEvid, file = "./plots/S6_weight_evidence_data.png",
       width = 8, height = 4, units = "in")


