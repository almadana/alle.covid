# Generates plots of the dynamics of the simulations. Generates a
# plot for the main text and a plot with a larger sample for the supplementary

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(segmented)
library(gridExtra)

set.seed(2691)
nPlotDyn <- 4

dynamicsData <- readRDS("./generated_data/SIR_dynamics_simulation.RDS") %>%
  as_tibble(.)
dynamicsFit <- readRDS("./generated_data/SIR_dynamics_fit.RDS") %>%
  as_tibble(.)

simulationsFit <- merge(dynamicsData, dynamicsFit, by=c("rep","allee"))
simulationsFit <- simulationsFit %>%
  mutate(., t_bp = ceiling(10^time.threshold),
         cumI_fit = ifelse(t<10^time.threshold,
                           10^(intercept + log10(t)*slopeI),
                           10^(intercept - log10(t_bp)*slopeF +log10(t_bp)*slopeI +log10(t)*slopeF))) %>%
  as_tibble(.)

##################################
#### scatter plot of initial slope and final slope for simulations
##################################
#plotChange <- dynamicsFit %>%
#  dplyr::mutate(., allee = c("w/o Allee", "Allee")[as.integer(allee)+1]) %>%
#  ggscatterhist(., x = "slopeI", y = "slopeF", color = "allee",
#                alpha = 0.5, size = 3, margin.plot = "boxplot",
#                palette = c("#b33018", "#14b74b"),
#                ggtheme = theme_bw(), xlab="Initial slope",
#                ylab="Slope after threshold",
#                legend.title = "",
#                margin.params = list(color = "allee"))


#############################
#### Plot the full dynamics of a selection of simulations
#############################
alleeSampleFit <- dplyr::filter(dynamicsFit, slopeRatio > 1 & allee)
slopeOrder <- base::order(alleeSampleFit$slopeRatio)
nSims <- length(slopeOrder)
ind <- round(c(seq(3, nSims, nSims/(nPlotDyn-1)), nSims-5))
alleeSamples <- alleeSampleFit$rep[ind]

nonAlleeSampleFit <- dplyr::filter(dynamicsFit, slopeRatio > 1 & !allee)
slopeOrder <- base::order(nonAlleeSampleFit$slopeRatio)
nSims <- length(slopeOrder)
ind <- round(c(seq(10, nSims, nSims/(nPlotDyn-1)), nSims-10))
nonAlleeSamples <- nonAlleeSampleFit$rep[ind]

plotsDynAllee <- dplyr::filter(simulationsFit,
                               allee == TRUE & rep %in% alleeSamples) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#b33018") +
  facet_wrap(~rep, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("with NPI") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="pink")

plotsDynNonAllee <- dplyr::filter(simulationsFit,
                                  allee == FALSE & rep %in% nonAlleeSamples) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#14b74b") +
  facet_wrap(~rep, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("without NPI") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="green")

dynsPlots <- grid.arrange(plotsDynAllee, plotsDynNonAllee, ncol=2,
                          left = "Cumulative infected")

ggsave("./plots/simDynamics.pdf", dynsPlots, width = 18, height = 11, units = "in")

##############################
######## make figure 1 #######
##############################

source("./1_Allee_phase_space.R")
fig1 <- ggarrange(allee1D, plotGrid, dynsPlots, nrow = 1, widths = c(2.8, 4, 3.7),
                  labels = "auto")
ggsave("./plots/fig1.png", fig1, width = 11, height = 3.6, units = "in")

################
# Plot extended sample of simulation dynamics for supplementary
################

suppN <- 40

sampleRepetitions <- sample(unique(dynamicsData$rep), suppN)
sampleSupp <- sampleRepetitions[order(sampleRepetitions)]

plotsDynAlleeSupp <- dplyr::filter(simulationsFit,
                               allee == TRUE & rep %in% sampleSupp) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#b33018") +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("with NPI") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="pink")


plotsDynNonAlleeSupp <- dplyr::filter(simulationsFit,
                                      allee == FALSE & rep %in% sampleSupp) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#14b74b") +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("without NPI") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="green")

dynsPlotsSupp <- grid.arrange(plotsDynAlleeSupp, plotsDynNonAlleeSupp, ncol=2,
                          left = "Cumulative infected")

ggsave("./plots/simDynamicsSupp.png", dynsPlotsSupp, width = 11, height = 15,
       units = "in")

