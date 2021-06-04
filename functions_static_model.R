library(tidyr)
library(gtable)

# function of effective social links by day
links_fun <- function(maxLinks, day){
  return(maxLinks * (day^4) / (day^4 + 7^4))
}

# calculate the expected number of links of detected individuals
expected_links <- function(maxLinks, pContact){
  days <- c(0:20)
  expectedLinks <- NULL
  for (n in c(1:length(pContact))) {
    pDay <- dgeom(days, pContact[n])
    pDay[length(pDay)] <- 1 - sum(pDay[1:(length(pDay)-1)])
    links <- links_fun(maxLinks, days)
    expectedLinks[n] <- sum(pDay * links)
  }
  return(expectedLinks)
}

# calculate Rt for given parameters
calculate_Rt <- function(maxDetected, I50, maxLinks, Infected, pCall,
                         NDailyCalls, pInfection, Npop){
  # Expected number of detected individuals
  pDetected <- maxDetected / (I50 + Infected)
  NDetected <- pDetected * Infected
  # Daily probability that a detected individual is contacted
  pContact <- 1 - (1-pCall)^(NDailyCalls/NDetected)
  # Number of expected social links for detected individuals
  detectedLinks <- expected_links(maxLinks, pContact)
  # Probability that a contact of an infected is susceptible
  pSusceptible<- (Npop - Infected - 1) / (Npop - 1)
  # Rt calculation
  Rt <- pSusceptible * pInfection * (pDetected*detectedLinks + (1-pDetected)*maxLinks)
  return(Rt)
}

# calculates Rt for a vector of values of Infected
calculate_Rt_range <- function(maxDetected, I50, maxLinks, Infected, pCall,
                               NDailyCalls, pInfection, Npop){
  RtVec <- vector()
  for (i in c(1:length(Infected))) {
    RtVec[i] <- calculate_Rt(maxDetected, I50, maxLinks, Infected[i], pCall,
                             NDailyCalls, pInfection, Npop)
  }
  return(RtVec)
}

##############################
# --- Functions that generate phase space -----
#############################

### Calculate Rt for a varying number of detection capacity
Rt_detected <- function(maxLinks, maxDetected, I50, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(maxDetected), nrow=length(infected))
  for(md in c(1:length(maxDetected))) {
    RtSpace[,md] <-  calculate_Rt_range(maxDetected=maxDetected[md],
                                        I50=I50[md],
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


### Calculate Rt for a varying number of daily calls
Rt_calls <- function(maxLinks, maxDetected, I50, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(NDailyCalls), nrow=length(infected))
  for(nc in c(1:length(NDailyCalls))) {
    RtSpace[,nc] <-  calculate_Rt_range(maxDetected=maxDetected, I50=I50,
                                        maxLinks=maxLinks, Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls[nc],
                                        pInfection=pInfection,
                                        Npop=Npop)
  }
  colnames(RtSpace) <- NDailyCalls
  rownames(RtSpace) <- infected
  return(RtSpace)
}


### Calculate Rt for a varying number of maximum contacts
Rt_links <- function(maxLinks, I50, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(maxLinks), nrow=length(infected))
  for(ml in c(1:length(maxLinks))) {
    RtSpace[,ml] <-  calculate_Rt_range(maxDetected=maxDetected,
                                        I50=I50,
                                        maxLinks=maxLinks[ml], Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls,
                                        pInfection=pInfection,
                                        Npop=Npop)
  }
  colnames(RtSpace) <- maxLinks
  rownames(RtSpace) <- infected
  return(RtSpace)
}

### Calculate Rt for a varying number of probability of infectious close contact
Rt_pInfection <- function(maxLinks, I50, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(pInfection), nrow=length(infected))
  for(pI in c(1:length(pInfection))) {
    RtSpace[,pI] <-  calculate_Rt_range(maxDetected=maxDetected,
                                        I50=I50,
                                        maxLinks=maxLinks, Infected=infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls,
                                        pInfection=pInfection[pI],
                                        Npop=Npop)
  }
  colnames(RtSpace) <- pInfection
  rownames(RtSpace) <- infected
  return(RtSpace)
}

# Enlongate and subsample matrix, and binarize, to draw raster in ggplot
# varNames is the names of the columns
# infectedSubsample indicates to keep 1 of every infectedSumsample infected values
enlongate_matrix <- function(inputMatrix, varNames, infectedSubsample) {
  longDF <- melt.array(inputMatrix, varnames = varNames) %>%
    dplyr::rename(., Rt = value) %>%
    dplyr::mutate(., Rt_exp = as.factor(as.integer(Rt > 1))) %>%
    as_tibble(.) %>%
    dplyr::filter(., (Infected %% infectedSubsample) == 0)
  return(longDF)
}

# plot 2D
plot_phase_space <- function(inputDF, reverse = FALSE){
  spacePlot <- ggplot(inputDF, aes(x = measure, y = Infected, fill = Rt_exp)) +
    geom_raster() +
    scale_fill_manual(values = c("#243faf", "#fa3d1b"),
                      name = element_blank(), labels = c(bquote(R[e]<1~"(Containment)"),
                                                         bquote(R[e]>1~" (Outbreak)"))) +
    scale_y_continuous(name = "Proportion infected", expand = c(0,0),
                       limits = c(0,.8), breaks = c(0, 0.4, 0.8)) +
    theme_bw()
  if (reverse) {
    spacePlot <- spacePlot +
      scale_x_reverse(name = xName, expand = c(0,0))
  }
  return(spacePlot)
}

