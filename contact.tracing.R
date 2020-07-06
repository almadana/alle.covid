###### Plot the state space of Rt against N infected for different interventions


# maxLinks - maximum number of close contacts that an
# infected can achieve 
# Infected - number of infected people in the population
# pInfection - probability that a close contact results 
# in an infectionof transmision through an effective
# link.

# pCall - probability that an attempt to reach a detected contact
# (e.g. phone call) succeeds 
# Npop is the total population size

library(tidyr)
library(ggplot2)
library(reshape)
library(dplyr)

# function of effective social links by day
links_fun <- function(maxLinks, day){
  return(maxLinks * (day^4) / (day^4 + 7^4))
}

# calculate the expected number of links of detected individuals
expected_links <- function(maxLinks, pContact){
  days <- c(1:20)
  pDay <- dgeom(days, pContact)
  links <- links_fun(maxLinks, days)
  expectedLinks <- sum(pDay * links)
  return(expectedLinks)
}

# calculate Rt for given parameters
calculate_Rt <- function(maxDetected, maxLinks, Infected, pCall, NDailyCalls,
                   pInfection, Npop){
  # Expected number of detected individuals
  pDetected <- maxDetected / (maxDetected + Infected)
  NDetected <- pDetected * Infected
  # Daily probability that a detected individual is contacted
  pContact <- 1 - (1-pCall)^(NDailyCalls/NDetected)
  # Number of expected social links for detected individuals
  detectedLinks <- expected_links(maxLinks, pContact)
  # Probability that an individual is susceptible
  pSusceptible<- (Npop - Infected) / Npop 
  # Rt calculation
  Rt <- pSusceptible * pInfection * (pDetected*detectedLinks + (1-pDetected)*maxLinks)
  return(Rt)
}

# calculates Rt for a vector of values of Infected
calculate_Rt_range <- function(maxDetected, maxLinks, Infected, pCall,
                               NDailyCalls, pInfection, Npop){
  RtVec <- vector()
  for (i in c(1:length(Infected))) {
    RtVec[i] <- calculate_Rt(maxDetected, maxLinks, Infected[i], pCall,
                             NDailyCalls, pInfection, Npop)
  }
  return(RtVec)
}

##############################
# Space Rt~ Infected + strategy

### Calculate Rt for a varying number of maximum contacts
Rt_calls <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(maxDetected), nrow=length(infected))
  for(md in c(1:length(maxDetected))) {
    RtSpace[,md] <-  calculate_Rt_range(maxDetected=maxDetected[md],
                                        maxLinks=maxLinks, Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls,
                                        pInfection=pInfection,
                                        Npop=Npop)
  }
  colnames(RtSpace) <- maxDetected
  rownames(RtSpace) <- infected
  return(RtSpace)
}

### Calculate Rt for a varying number of maximum links
Rt_links <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop){
  RtLinks <- matrix(NA, ncol=length(maxLinks), nrow=Npop)
  infected <- c(1:Npop)
  for(ml in c(1:length(maxLinks))) {
    RtLinks[,ml] <-  calculate_Rt_range(maxDetected=maxDetected,
                                        maxLinks=maxLinks[ml], Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls,
                                        pInfection=pInfection,
                                        Npop=Npop)
  }
  return(RtLinks)
}

### Calculate Rt for a varying number of pInfection
Rt_pInfection <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop){
  RtInf <- matrix(NA, ncol=length(pInfection), nrow=Npop)
  infected <- c(1:Npop)
  for(pI in c(1:length(pInfection))) {
    RtInf[,pI] <-  calculate_Rt_range(maxDetected=maxDetected,
                                        maxLinks=maxLinks, Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls,
                                        pInfection=pInfection[pI],
                                        Npop=Npop)
  }
  return(RtInf)
}


###############################
# Make the space plots
###############################

#general parameters
nTreeDay[nTreeDay < 1] <- 1
pCall <- 0.2
pInfection <- 0.3
Npop <- 2000
maxLinks <- 12
propInfectedMax = 0.8

# Space plot for maximum number of people that can be detected
nByTree <- 10
nTreeDay <- c(1:700) * 0.2
maxDetected_Vec <- nByTree * nTreeDay
callsMatrix <- Rt_calls(maxLinks = maxLinks, maxDetected = maxDetected_Vec,
                        pCall = pCall, pInfection = pInfection,
                        NDailyCalls = NDailyCalls, Npop = Npop,
                        propInfectedMax = propInfectedMax)

# reshape matrix into dataframe for plotting
callsLong <- melt.array(callsMatrix, varnames = c("Infected", "maxDetected")) %>%
  dplyr::rename(., Rt = value) %>%
  dplyr::mutate(., Rt_exp = as.factor(as.integer(Rt > 1))) %>%
  as_tibble(.) %>%
  dplyr::filter(., (Infected %% 2) == 0)

RtSpacePlot <- ggplot(callsLong, aes(x = maxDetected, y = Infected, fill = Rt_exp)) +
  geom_raster() +
  scale_fill_manual(values = c("#243faf", "#fa3d1b"),
                    name = "Rt", labels = c("< 1", ">= 1")) +
  scale_x_continuous(name = "Detection capacity (N infected)", expand = c(0,0)) +
  scale_y_continuous(name = "Infected individuals", expand = c(0,0)) +
#  geom_segment(data = arrowDf,
#               aes(x = x, xend = x, y = y, yend = y+arrowLength, fill = "black"),
#            arrow = arrow(length = unit(0.5, "cm")))+
  annotate("text", x = 500, y = 750, label = "Growdth", size = 20) +
  annotate("text", x = 950, y = 200, label = "Containment", size = 20,
           color = "white") +
  theme_bw()


