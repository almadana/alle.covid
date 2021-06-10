###### Plot the needed number of vaccinated people to contain
# epidemic for different levels of NPIs

# Variables
#
# Npop: overall population
# pCall: probability that an attempt to contact individual is succesful
# pInfection: probability that a close contact is infected
# maxLinks: average number of close contacts that non detected infected individuals make
# NDailyCalls: daily number of contacting attempts that can be performed
# maxDetected: maximum number of infected individuals that can be detected
# infectedSubsample: plotting parameter to subsample infected for 2D rasters
# propInfectedMax: plotting parameter to truncate 2D plots vertically

library(tidyr)
library(ggplot2)
library(reshape)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)
library(gridExtra)
library(ggpubr)

source("./functions_static_model.R")

###############################
# ---- Make the space plots -----
###############################


# Function to compute HI threshold for different strength of NPI
find_HI_threshold <- function(longNPITable, NPI_col, Npop) {
  thresholds <- longNPITable %>%
    dplyr::mutate(., Psuc1=(Npop-Infected-1)/(Npop-1),
                  undo_suc_Rt=Rt/Psuc1,
                  vacSoRt1=(Npop-Infected-1-(Npop-1)/undo_suc_Rt)/Npop) %>%
  dplyr::group_by_at(., NPI_col) %>%
  dplyr::summarize(., maxVac=max(vacSoRt1)) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(., HI_threshold=max(maxVac,0)) %>%
  dplyr::ungroup(.)
  return(thresholds)
}

# general parameters
pCall <- 0.1
pInfection <- 0.2
Npop <- 2000
maxLinks <- 14
propInfectedMax <- 1
NDailyCalls <- 800
infectedSubsample <- 3
maxDetected <- 200

# ranges for 
maxDetected_Vec <- seq(0, 1000, 5)
dailyCalls_Vec <- seq(20, 1500, 5)
maxLinks_Vec <- seq(2, 30, 0.1)
pInfection_Vec <- seq(0.1, 0.3, 0.0015)


## Plot HI threshold for the strenghts of individual NPIs
#detectedHI <- Rt_detected(maxLinks=maxLinks, maxDetected=maxDetected_Vec,
#                            I50=maxDetected_Vec, pCall=pCall,
#                            pInfection=pInfection,
#                        NDailyCalls=NDailyCalls, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>%
#                enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
#                find_HI_threshold(., "maxDetected", Npop=Npop)
#
#callsHI <- Rt_calls(maxLinks=maxLinks, maxDetected=maxDetected,
#                      I50=maxDetected, pCall=pCall,
#                      pInfection=pInfection,
#                        NDailyCalls=dailyCalls_Vec, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>%
#            enlongate_matrix(., c("Infected", "maxCalls"), infectedSubsample) %>%
#            find_HI_threshold(., "maxCalls", Npop=Npop)
#
#linksHI <- Rt_links(maxLinks=maxLinks_Vec, maxDetected=maxDetected,
#                        I50=maxDetected, pCall=pCall, pInfection=pInfection,
#                        NDailyCalls=NDailyCalls, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>% 
#            enlongate_matrix(., c("Infected", "maxLinks"), infectedSubsample) %>%
#            find_HI_threshold(., "maxLinks", Npop=Npop)
#
#infectionHI <- Rt_pInfection(maxLinks=maxLinks, maxDetected=maxDetected,
#                               I50=maxDetected, pCall=pCall,
#                               pInfection=pInfection_Vec,
#                        NDailyCalls=NDailyCalls, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>%
#            enlongate_matrix(., c("Infected", "pInfection"), infectedSubsample) %>%
#            find_HI_threshold(., "pInfection", Npop=Npop)
#
#
#detectedHIPlot <- ggplot(detectedHI, aes(x=maxDetected, y=HI_threshold*100)) +
#  geom_line(size=1) +
#  xlab("Maximum detection capacity (K)") +
#  ylab("Vaccination HI threshold (%)") +
#  theme_classic()
#
#callsHIPlot <- ggplot(callsHI, aes(x=maxCalls, y=HI_threshold*100)) +
#  geom_line(size=1) +
#  xlab("Maximum call capacity") +
#  ylab("Vaccination HI threshold (%)") +
#  theme_classic()
#
#linksHIPlot <- ggplot(linksHI, aes(x=maxLinks, y=HI_threshold*100)) +
#  geom_line(size=1) +
#  xlab("Maximum social links (Lmax)") +
#  ylab("Vaccination HI threshold (%)") +
#  theme_classic()
#
#infectionHIPlot <- ggplot(infectionHI, aes(x=pInfection, y=HI_threshold*100)) +
#  geom_line(size=1) +
#  xlab("Link infection probability (Plink)") +
#  ylab("Vaccination HI threshold (%)") +
#  theme_classic()
#
## Put plots for different NPIs together
#plotlist = list(detectedHIPlot, callsHIPlot, linksHIPlot, infectionHIPlot)
#allMeasuresPlots <- ggpubr::ggarrange(plotlist=plotlist, ncol=2, nrow=2)
#
#ggsave("./plots/vaccination_NPI1.png", allMeasuresPlots, width=15,
#       height=15, units="cm")
#

### Plot for detection capacity the HI threshold when changing other
# parameters

# general parameters
pCall <- 0.1
pInfection <- 0.2
Npop <- 2000
maxLinks <- 14
propInfectedMax <- 1
NDailyCalls <- 800
infectedSubsample <- 3

maxDetected_Vec <- seq(0, 1000, 5)

## ranges of daily calls
#dailyCalls_Vec <- seq(300, 1500, 500)
#detected_Call_HI <- NULL
#for (dc in c(1:length(dailyCalls_Vec))) {
#  tempHIDf <- Rt_detected(maxLinks=maxLinks, maxDetected=maxDetected_Vec,
#                              I50=maxDetected_Vec, pCall=pCall,
#                              pInfection=pInfection,
#                          NDailyCalls=dailyCalls_Vec[dc], Npop=Npop,
#                          propInfectedMax=propInfectedMax) %>%
#                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
#                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
#                  dplyr::mutate(., nCalls=dailyCalls_Vec[dc])
#  detected_Call_HI <- rbind(detected_Call_HI, tempHIDf)
#}
#detected_Call_HI$nCalls_F <- factor(detected_Call_HI$nCalls)
#
#detected_Call_Plot <- detected_Call_HI %>%
#  dplyr::filter(., HI_threshold>0) %>%
#  ggplot(., aes(x=maxDetected, y=HI_threshold*100, color=nCalls_F)) +
#  geom_line(size=1) +
#  xlab("Maximum detection capacity (K)") +
#  ylab("Vaccination HI threshold (%)") +
#  labs(color="Max calls") +
#  theme_classic()


# ranges of max links
maxLinks_Vec <- seq(10, 25, 5)
detected_Link_HI <- NULL
for (ml in c(1:length(maxLinks_Vec))) {
  tempHIDf <- Rt_detected(maxLinks=maxLinks_Vec[ml],
                          maxDetected=maxDetected_Vec,
                          I50=maxDetected_Vec, pCall=pCall,
                          pInfection=pInfection,
                          NDailyCalls=NDailyCalls, Npop=Npop,
                          propInfectedMax=propInfectedMax) %>%
                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
                  dplyr::mutate(., maxLinks=maxLinks_Vec[ml])
  detected_Link_HI <- rbind(detected_Link_HI, tempHIDf)
}
detected_Link_HI$maxLinks_F <- factor(detected_Link_HI$maxLinks)

detected_Link_Plot <- detected_Link_HI %>%
  dplyr::filter(., HI_threshold>0) %>%
  ggplot(., aes(x=maxDetected, y=HI_threshold*100, color=maxLinks_F)) +
  geom_line(size=0.6) +
  xlab("Maximum detection capacity (K)") +
  ylab("Vaccination HI threshold (%)") +
  labs(color="Max links") +
  theme_classic()


# Ranges of pInfection
pInfection_Vec <- seq(0.15, 0.3, 0.05)
detected_Pinfection_HI <- NULL
for (pinf in c(1:length(pInfection_Vec))) {
  tempHIDf <- Rt_detected(maxLinks=maxLinks,
                          maxDetected=maxDetected_Vec,
                          I50=maxDetected_Vec, pCall=pCall,
                          pInfection=pInfection_Vec[pinf],
                          NDailyCalls=NDailyCalls, Npop=Npop,
                          propInfectedMax=propInfectedMax) %>%
                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
                  dplyr::mutate(., pInfection=pInfection_Vec[pinf])
  detected_Pinfection_HI <- rbind(detected_Pinfection_HI, tempHIDf)
}
detected_Pinfection_HI$pInfection_F <- factor(detected_Pinfection_HI$pInfection)

detected_Pinfection_Plot <- detected_Pinfection_HI %>%
  dplyr::filter(., HI_threshold>0) %>%
  ggplot(., aes(x=maxDetected, y=HI_threshold*100, color=pInfection_F)) +
  geom_line(size=0.6) +
  xlab("Maximum detection capacity (K)") +
  ylab("Vaccination HI threshold (%)") +
  labs(color=bquote(~b[link])) +
  theme_classic()

plotlist2 = list(detected_Link_Plot, detected_Pinfection_Plot)
interactionMeasuresPlot <- ggpubr::ggarrange(plotlist=plotlist2, ncol=2, nrow=1)

ggsave("./plots/vaccination_NPI2.png", interactionMeasuresPlot, width=20,
       height=10, units="cm")


######## Plot 1D phase space for differente vaccination proportions

# general parameters
pCall <- 0.1
pInfection <- 0.2
Npop <- 2000
maxLinks <- 14
propInfectedMax <- 1
NDailyCalls <- 800
infectedSubsample <- 3
maxDetected <- 200
Infected <- c(1:2000)

vacProp <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
maxDetected_logistic <- 0
I50_w <- 400

vacNPIdf <- data.frame()
for (vp in c(1:length(vacProp))) {
  logisticR <- calculate_Rt(maxDetected=maxDetected_logistic, I50=maxDetected,
                        maxLinks=maxLinks, Infected=Infected,
                        pCall=pCall, NDailyCalls=NDailyCalls,
                        pInfection=pInfection, Npop=Npop,
                        propVac=vacProp[vp])
  logisticDf <- data.frame(Infected=Infected, Rt=logisticR,
                           type="No-NPI", Vaccinated=vacProp[vp]*100)
  weakAlleeR <- calculate_Rt(maxDetected=maxDetected, I50=I50_w,
                        maxLinks=maxLinks, Infected=Infected,
                        pCall=pCall, NDailyCalls=NDailyCalls,
                        pInfection=pInfection, Npop=Npop,
                        propVac=vacProp[vp])
  weakAlleeDf <- data.frame(Infected=Infected, Rt=weakAlleeR,
                           type="Weak-Allee", Vaccinated=vacProp[vp]*100)
  alleeR <- calculate_Rt(maxDetected=maxDetected, I50=maxDetected,
                        maxLinks=maxLinks, Infected=Infected,
                        pCall=pCall, NDailyCalls=NDailyCalls,
                        pInfection=pInfection, Npop=Npop,
                        propVac=vacProp[vp])
  alleeDf <- data.frame(Infected=Infected, Rt=alleeR,
                           type="Strong-Allee", Vaccinated=vacProp[vp]*100)
  vacNPIdf <- rbind(vacNPIdf, logisticDf, weakAlleeDf, alleeDf)
}


weakAlleeVacPlot <- dplyr::filter(vacNPIdf, type=="Weak-Allee" & Rt>0) %>%
  dplyr::mutate(., Vaccinated=factor(Vaccinated)) %>%
  ggplot(., aes(x=Infected, y=Rt, color=Vaccinated)) +
  geom_line(size = 1)

logisticAlleeVacPlot <- dplyr::filter(vacNPIdf, type=="No-NPI" & Rt>0) %>%
  dplyr::mutate(., Vaccinated=factor(Vaccinated)) %>%
  ggplot(., aes(x=Infected, y=Rt, color=Vaccinated)) +
  geom_line(size = 1)

alleeVacPlot <- dplyr::filter(vacNPIdf, type=="Strong-Allee" & Rt>0) %>%
  dplyr::mutate(., Vaccinated=factor(Vaccinated)) %>%
  ggplot(., aes(x=Infected, y=Rt, color=Vaccinated)) +
  geom_line(size=0.6) +
  geom_hline(yintercept=1, size=1, linetype="dashed") +
  xlab("Number infected") +
  ylab(bquote('Reproductive number ' ~R[e])) +
#  scale_x_continuous(breaks=c(0, 1000, 2000), labels=c(0, 0.5, 1)) +
  theme_classic()


fig2 <- ggarrange(interactionMeasuresPlot, alleeVacPlot, nrow=1,labels="AUTO",
                  widths=c(0.65, 0.35))
ggsave("./plots/fig2.png", fig2, width=10, height=3, units="in")

