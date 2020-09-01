library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggrepel)
library(class)

nExamples <- 8

# load data
countiesFit <- readRDS("./generated_data/countiesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI)

simFit <- readRDS("./generated_data/SIR_dynamics_fit_parameter_exploration.RDS") %>%
  dplyr::mutate(allee = factor(ifelse(allee,"with NPI","without NPI")))


countiesFit$allee = factor("U.S. counties", levels = c("U.S. counties", levels(simFit$allee)))
simFit$Id <- seq(990000, 990000 + nrow(simFit)-1, 1)
 
mergeCols <- c("Id", "allee", "slopeI", "slopeF", "slopeRatio")
#("Id","allee","slopeI","slopeF","slopeRatio")
plotTitle <- paste("p=", sim_p, ", I50=", sim_I50, ", betaMax=", sim_betaMax,
                   ", gammaMax=", sim_gammaMax,
                   ", dispersion=", sim_dispersion, sep = "")
annotationDf <- data.frame(x = rep(0.31, 5),
                           y = seq(from=2, to=1.3, length.out=5))

plot_slopes <- list()
proportionSimilarAllee <- NULL
for (sim in unique(simFit$simN)) {
  parameterFit <- dplyr::filter(simFit, simN == sim)
  sim_p <- paste("p ==", as.character(parameterFit$p.x[1]))
  sim_I50 <- paste("I[50] ==", as.character(parameterFit$I50[1]))
  sim_dispersion <- paste("dispersion ==", as.character(parameterFit$dispersion[1]))
  sim_betaMax <- paste("beta[Max] ==", as.character(parameterFit$betaMax[1]))
  sim_gammaMax <- paste("gamma[Max] ==", as.character(parameterFit$gammaMax[1]))
  parameterDf <- annotationDf
  parameterDf$parText <- c(sim_p, sim_I50, sim_dispersion, sim_betaMax, sim_gammaMax)

  sim_countiesFit <- merge(parameterFit[,mergeCols], countiesFit, by=mergeCols, all=T) %>%
    as_tibble(.)

  sim_countiesFit$breakout <- sim_countiesFit$slopeRatio > 1
  sim_countiesFit$breakout_allee = interaction(sim_countiesFit$allee, sim_countiesFit$breakout)
  sim_countiesFit$allee <- as.character(sim_countiesFit$allee)
  sim_countiesFit$allee <- factor(sim_countiesFit$allee,
                                  levels = c("U.S. Counties", "with NPI", "without NPI"))

  xlims=c(0,.39)
  ylims= c(-.01,2)

  plot_slopes[[sim]] <- sim_countiesFit %>%
    #filter(slopeRatio>1) %>%
    ggplot(., aes(group=allee, fill=allee, x=log10(1+slopeI), y=log10(1+slopeRatio))) +
    stat_density_2d(geom="polygon", aes(fill=allee, alpha=..level..), contour=T, bins=14) + 
    geom_point(shape=20, data=(countiesFit), alpha=1,
               aes(col=allee), size = 0.8) +
    scale_fill_manual(values = c("#FFFFFF","#b33018","#14b74b")) +
    theme_classic() +
    scale_color_manual(values=c("#000080","#b33018","#14b74b"))+
    geom_abline(slope=0,intercept = log10(2),linetype="dotted") +
#    geom_text_repel(data=sim_countiesFit, aes(label=label), size=3) +
    #annotate("text", x=0.03, y=log10(2)+0.1, label="equal slopes", size=3) +
#    geom_text(data = parameterDf, aes(label=parText, x = x, y = y), size=3) +
    xlim(xlims) +
    ylim(ylims) +
  #  labs(x="log10(1+slope before breakout point)",y="log10(1+slope after/slope before)")
    labs(x=expression(paste(Log10, "(", 1+mu[1], ")")),
         y=expression(paste(Log10, "(", 1+frac(mu[2],mu[1]), ")"))) +
#    ggtitle(plotTitle) +
    theme(legend.position = "top") +  
    guides(alpha="none",shape="none",colour="none", fill="none")

    for (n in c(1:nrow(parameterDf))) {
      plot_slopes[[sim]] <- plot_slopes[[sim]] +
        annotate("text", label = parameterDf$parText[n], x = parameterDf$x[n],
                 y = parameterDf$y[n], size = 3.5, parse = TRUE)
    }

  plot_slopes[[sim]]$labels$fill=""

  # use KNN to see proportion of counties with closest neighbor Allee
  trainPredictors <- dplyr::select(parameterFit, slopeI, slopeRatio)
  testPredictors <- dplyr::select(countiesFit, slopeI, slopeRatio)
  trainPredictors <- dplyr::mutate(trainPredictors, slopeI = log10(1+slopeI),
                                   slopeRatio = log10(1+slopeRatio))
  testPredictors <- dplyr::mutate(testPredictors, slopeI = log10(1+slopeI),
                                   slopeRatio = log10(1+slopeRatio))
  intLabel <- as.integer(parameterFit$allee == "with NPI")
  trainLabels <- c("non Allee", "Allee")[intLabel + 1]
  naInds <- which(is.na(trainPredictors$slopeRatio))
  infInds <- which(is.infinite(trainPredictors$slopeRatio))
  removeInds <- c(naInds, infInds)
  if (length(removeInds)>0) {
    trainPredictors <- trainPredictors[-removeInds,]
    trainLabels <- trainLabels[-removeInds]
  }
  naIndsTest <- which(is.na(testPredictors$slopeRatio))
  if (length(naIndsTest)>0) {
    testPredictors <- testPredictors[-naIndsTest,]
  }
  knnResult <- knn(train = trainPredictors, test = testPredictors,
                   cl = trainLabels, k = 10)
  proportionSimilarAllee[length(proportionSimilarAllee)+1] <-
    mean(knnResult == "Allee")
}

#plot density plots
nSplits <- 4
plotInds <- c(1:72)
nCol <- 3
chunkSize <- max(plotInds)/nSplits
indSplits <- split(plotInds, ceiling(seq_along(plotInds)/chunkSize))

figWidth <- 10
figHeight <- 15
for (n in 1:length(indSplits)) {
  nRow <- ceiling(length(indSplits[[n]])/nCol)
  split <- ggarrange(plotlist = plot_slopes[indSplits[[n]]],
                     ncol = nCol, nrow = nRow)
  fileName <- paste("./plots/parameter_exploration",
                    as.character(n), ".png", sep="")
  ggsave(split, file=fileName, width=figWidth,
         height=figHeight, units = "in")
}

# plot histogram with similarities
proportionSimilarAllee <- data.frame(Allee_class = proportionSimilarAllee)
knnPlot <- ggplot(data = proportionSimilarAllee, aes(x = Allee_class)) +
    geom_histogram(color = "gray", fill="#b33018") +
    xlab("Proportion of dynamics classified as Allee") +
    ylab("Parameter combinations") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5, linetype="solid"),
          axis.line.y = element_line(size=0.5, linetype="solid")) +
    xlim(0, 1)

ggsave(knnPlot, file="./plots/parameter_exploration_KNN_histogra.png",
       width=4, height=4, units = "in")


# Plot proportion of simulations with weight of evidence > 0.95
proportionSegmented <- group_by(simFit, simN, allee) %>%
  summarize(., highlySegmented = mean(weighted.evidence > 0.95))

plotSegmentedProportion <- proportionSegmented %>%
  ggplot(., aes(highlySegmented, fill = allee)) +
  geom_histogram(color = "gray") +
  facet_wrap(~allee, ncol=1) +
  xlab("Proportion with weight of evidence > 0.95") +
  ylab("Parameter combinations") +
  scale_fill_manual(values = c("#b33018","#14b74b")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid"))


ggsave(plotSegmentedProportion, file="./plots/Supp_proportionSegmented.png",
       width = 4, height = 4, units = "in")

