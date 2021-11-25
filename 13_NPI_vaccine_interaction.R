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
library(viridis)

source("./functions_static_model.R")

###############################
# ---- Make the space plots -----
###############################

# Function to compute HI threshold for different strength of NPI
find_HI_threshold <- function(longNPITable, NPI_col, Npop) {
  # First, we calculate X = blink(L*fq+Lmax*fnq)
  # Then, using Psuc = S/(Npop-1), Pnvs=(Nnv-I)/S, Pnv=Nnv/(Npop-1)
  # Pv = 1-Pnv
  # where nv = non vaccinated, and Pnvs is the probability
  # of a susceptible not being vaccinated. Note that we
  # assume that all the Infected were not vaccinated, so
  # Nnv-I is the number of non vaccinated susceptibles.
  # Then, Re = Psuc*Pnvs*X
  # Solving for Re=1, we have
  # 1=[S/(Npop-1)]*[(Nnv-I)/S]*X =
  # (Nnv-I)/(Npop-1)*X = (Pnv-(1-Psuc))*X = [(1-Pv)-(1-Psuc)]*X =
  # (Psuc-Pv)*X  ---->   Pv = Psuc - 1/X
  thresholds <- longNPITable %>%
    dplyr::mutate(., Psuc1=(Npop-Infected-1)/(Npop-1),
                  undoSucRt=Rt/Psuc1, # same as blink(L*fq+Lmax*fnq), or X
                  vacSoRt1=Psuc1-1/undoSucRt) %>%
                  #vacSoRt1=(Npop-Infected-1-(Npop-1)/undoSucRt)/Npop) %>%
                  #vacSoRt1=1-1/Rt) %>%
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
I50 <- 200

# ranges for 
maxDetectedVec <- seq(0, I50, 5)
maxLinksVec <- seq(10, 25, 5)
pInfectionVec <- seq(0.15, 0.3, 0.05)

### Plot for detection capacity the HI threshold when changing other
# parameters

# ranges of max links
detectedLink_HI <- NULL
for (ml in c(1:length(maxLinksVec))) {
  tempHIDf <- Rt_detected(maxLinks=maxLinksVec[ml],
                          maxDetected=maxDetectedVec,
                          I50=I50, pCall=pCall,
                          pInfection=pInfection,
                          NDailyCalls=NDailyCalls, Npop=Npop,
                          propInfectedMax=propInfectedMax) %>%
                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
                  dplyr::mutate(., maxLinks=maxLinksVec[ml])
  detectedLink_HI <- rbind(detectedLink_HI, tempHIDf)
}
detectedLink_HI$maxLinks_F <- factor(detectedLink_HI$maxLinks)

detectedLinkPlot <- detectedLink_HI %>%
  dplyr::filter(., HI_threshold>0) %>%
  ggplot(., aes(x=maxDetected, y=HI_threshold*100, group=maxLinks, color=maxLinks_F)) +
  geom_line(size=0.8) +
  scale_color_viridis(discrete=TRUE) +
  xlab("Maximum detection capacity (K)") +
  ylab("Vaccination HI threshold (%)") +
  labs(color="Max links") +
  theme_classic()


# Ranges of pInfection
pInfectionVec <- seq(0.15, 0.3, 0.05)
detected_Pinfection_HI <- NULL
for (pinf in c(1:length(pInfectionVec))) {
  tempHIDf <- Rt_detected(maxLinks=maxLinks,
                          maxDetected=maxDetectedVec,
                          I50=I50, pCall=pCall,
                          pInfection=pInfectionVec[pinf],
                          NDailyCalls=NDailyCalls, Npop=Npop,
                          propInfectedMax=propInfectedMax) %>%
                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
                  dplyr::mutate(., pInfection=pInfectionVec[pinf])
  detected_Pinfection_HI <- rbind(detected_Pinfection_HI, tempHIDf)
}
detected_Pinfection_HI$pInfection_F <- factor(detected_Pinfection_HI$pInfection)

detected_PinfectionPlot <- detected_Pinfection_HI %>%
  dplyr::filter(., HI_threshold>0) %>%
  ggplot(., aes(x=maxDetected, y=HI_threshold*100, color=pInfection_F)) +
  geom_line(size=0.8) +
  scale_color_viridis(discrete=TRUE) +
  xlab("Maximum detection capacity (K)") +
  ylab("Vaccination HI threshold (%)") +
  labs(color=bquote(~b[link])) +
  theme_classic()

plotlist2 = list(detectedLinkPlot, detected_PinfectionPlot)
interactionMeasuresPlot <- ggpubr::ggarrange(plotlist=plotlist2, ncol=2, nrow=1)

ggsave("./plots/vaccination_NPI2.pdf", interactionMeasuresPlot, width=20,
       height=10, units="cm")



## ranges of daily calls
#dailyCallsVec <- seq(300, 1500, 500)
#detected_Call_HI <- NULL
#for (dc in c(1:length(dailyCallsVec))) {
#  tempHIDf <- Rt_detected(maxLinks=maxLinks, maxDetected=maxDetectedVec,
#                              I50=maxDetectedVec, pCall=pCall,
#                              pInfection=pInfection,
#                          NDailyCalls=dailyCallsVec[dc], Npop=Npop,
#                          propInfectedMax=propInfectedMax) %>%
#                  enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
#                  find_HI_threshold(., "maxDetected", Npop=Npop) %>%
#                  dplyr::mutate(., nCalls=dailyCallsVec[dc])
#  detected_Call_HI <- rbind(detected_Call_HI, tempHIDf)
#}
#detected_Call_HI$nCalls_F <- factor(detected_Call_HI$nCalls)
#
#detected_CallPlot <- detected_Call_HI %>%
#  dplyr::filter(., HI_threshold>0) %>%
#  ggplot(., aes(x=maxDetected, y=HI_threshold*100, color=nCalls_F)) +
#  geom_line(size=1) +
#  xlab("Maximum detection capacity (K)") +
#  ylab("Vaccination HI threshold (%)") +
#  labs(color="Max calls") +
#  theme_classic()

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
  ggplot(., aes(x=Infected/Npop, y=Rt, color=Vaccinated)) +
  geom_line(size=0.6) +
  geom_hline(yintercept=1, size=1, linetype="dashed") +
  geom_line(size=0.8) +
  scale_color_viridis(discrete=TRUE, direction=-1) +
  xlab("Proportion infected") +
  ylab(bquote('Reproductive number ' ~R[e])) +
#  scale_x_continuous(breaks=c(0, 1000, 2000), labels=c(0, 0.5, 1)) +
  theme_classic()


fig2 <- ggarrange(interactionMeasuresPlot, alleeVacPlot, nrow=1,labels="AUTO",
                  widths=c(0.65, 0.35))
ggsave("./plots/fig2.pdf", fig2, width=10, height=3, units="in")


## Plot HI threshold for the strenghts of individual NPIs
#detectedHI <- Rt_detected(maxLinks=maxLinks, maxDetected=maxDetectedVec,
#                            I50=maxDetectedVec, pCall=pCall,
#                            pInfection=pInfection,
#                        NDailyCalls=NDailyCalls, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>%
#                enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample) %>%
#                find_HI_threshold(., "maxDetected", Npop=Npop)
#
#callsHI <- Rt_calls(maxLinks=maxLinks, maxDetected=maxDetected,
#                      I50=maxDetected, pCall=pCall,
#                      pInfection=pInfection,
#                        NDailyCalls=dailyCallsVec, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>%
#            enlongate_matrix(., c("Infected", "maxCalls"), infectedSubsample) %>%
#            find_HI_threshold(., "maxCalls", Npop=Npop)
#
#linksHI <- Rt_links(maxLinks=maxLinksVec, maxDetected=maxDetected,
#                        I50=maxDetected, pCall=pCall, pInfection=pInfection,
#                        NDailyCalls=NDailyCalls, Npop=Npop,
#                        propInfectedMax=propInfectedMax) %>% 
#            enlongate_matrix(., c("Infected", "maxLinks"), infectedSubsample) %>%
#            find_HI_threshold(., "maxLinks", Npop=Npop)
#
#infectionHI <- Rt_pInfection(maxLinks=maxLinks, maxDetected=maxDetected,
#                               I50=maxDetected, pCall=pCall,
#                               pInfection=pInfectionVec,
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

