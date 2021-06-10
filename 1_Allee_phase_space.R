###### Plot the state space of Rt against N infected for different interventions

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

#########
# ---- Basic model functions ----
########

# function of effective social links by day
links_fun <- function(maxLinks, day){
  return(maxLinks * (day^4) / (day^4 + 7^4))
}

# calculate the expected number of links of detected individuals
expected_links <- function(maxLinks, pContact){
  days <- c(0:20)
  pDay <- dgeom(days, pContact)
  pDay[length(pDay)] <- 1 - sum(pDay[1:(length(pDay)-1)])
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
  # Probability that a contact of an infected is susceptible
  pSusceptible<- (Npop - Infected - 1) / (Npop - 1)
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
# --- Functions that generate phase space -----
#############################

### Calculate Rt for a varying number of detection capacity
Rt_detected <- function(maxLinks, maxDetected, pCall,  pInfection,
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


### Calculate Rt for a varying number of daily calls
Rt_calls <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(NDailyCalls), nrow=length(infected))
  for(nc in c(1:length(NDailyCalls))) {
    RtSpace[,nc] <-  calculate_Rt_range(maxDetected=maxDetected,
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
Rt_links <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(maxLinks), nrow=length(infected))
  for(ml in c(1:length(maxLinks))) {
    RtSpace[,ml] <-  calculate_Rt_range(maxDetected=maxDetected,
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
Rt_pInfection <- function(maxLinks, maxDetected, pCall,  pInfection,
                     NDailyCalls, Npop, propInfectedMax = 1){
  infected <- c(1:round(Npop*propInfectedMax))
  RtSpace <- matrix(NA, ncol=length(pInfection), nrow=length(infected))
  for(pI in c(1:length(pInfection))) {
    RtSpace[,pI] <-  calculate_Rt_range(maxDetected=maxDetected,
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


###############################
# ---- Make the space plots -----
###############################

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

# Space plot for maximum number of people that can be detected
detectedLong <- Rt_detected(maxLinks = maxLinks, maxDetected = maxDetected_Vec,
                        pCall = pCall, pInfection = pInfection,
                        NDailyCalls = NDailyCalls, Npop = Npop,
                        propInfectedMax = propInfectedMax) %>%
                enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample)

detectedPlot <- dplyr::rename(detectedLong, measure = maxDetected) %>%
  dplyr::mutate(Infected =Infected/ Npop) %>% 
  plot_phase_space(.) +
  ggtitle("Tracing capacity") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.title.x = element_text(size = 10)) +
  scale_x_continuous(breaks = c(0, 400, 800),
                     name = bquote("Detection capacity (K)"),
                     expand = c(0,0))

# Space plot for maximum number of calls (speed of detection)
callsLong <- Rt_calls(maxLinks = maxLinks, maxDetected = maxDetected,
                        pCall = pCall, pInfection = pInfection,
                        NDailyCalls = dailyCalls_Vec, Npop = Npop,
                        propInfectedMax = propInfectedMax) %>%
            enlongate_matrix(., c("Infected", "maxCalls"), infectedSubsample)

callsPlot <- dplyr::rename(callsLong, measure = maxCalls) %>%
  dplyr::mutate(Infected =Infected/ Npop) %>% 
  plot_phase_space(.) +
  ggtitle("Tracing speed") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  scale_x_continuous(breaks = c(0, 500, 1000),
                     name = bquote("Maximum daily calls ("*N[calls]*")"),
                     expand = c(0,0))

# Space plot for maximum number of links
linksLong <- Rt_links(maxLinks = maxLinks_Vec, maxDetected = maxDetected,
                        pCall = pCall, pInfection = pInfection,
                        NDailyCalls = NDailyCalls, Npop = Npop,
                        propInfectedMax = propInfectedMax) %>% 
            enlongate_matrix(., c("Infected", "maxLinks"), infectedSubsample)
linksPlot <- dplyr::rename(linksLong, measure = maxLinks) %>%
  dplyr::mutate(Infected =Infected/ Npop) %>% 
  plot_phase_space(.) +
  ggtitle("Social distancing") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  scale_x_continuous(breaks = c(10, 20),
                     name = bquote("Maximum social links ("*L[max]*")"),
                     expand = c(0,0))

# Space plot for pInfection 
infectionLong <- Rt_pInfection(maxLinks = maxLinks, maxDetected = maxDetected,
                        pCall = pCall, pInfection = pInfection_Vec,
                        NDailyCalls = NDailyCalls, Npop = Npop,
                        propInfectedMax = propInfectedMax) %>%
            enlongate_matrix(., c("Infected", "pInfection"), infectedSubsample)
infectionPlot <- dplyr::rename(infectionLong, measure = pInfection) %>%mutate(Infected =Infected/ Npop) %>% 
  plot_phase_space(.) +
  ggtitle("Hygiene measures") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) +
  scale_x_continuous(breaks = c(0.15, 0.25),
                     name = bquote("Contact infectiousness ("*b[link]*")"),
                     expand = c(0,0))



# ---- arrange plots into grid -----
legendPlot <- detectedPlot + theme(legend.position = "bottom",
                                   legend.key.size = unit(.5,"cm"),
                                   legend.text = element_text(size = 9))
legend <- gtable_filter(ggplotGrob(legendPlot), "guide-box")

space2Dplots <- grid.arrange(
                     arrangeGrob(detectedPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 callsPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 infectionPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 linksPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 ncol=2, left = "Proportion infected")#,legend,
                     #heights=c(10, .5)
                     )


ggsave("./plots/phase_space.png", space2Dplots, width = 15, height = 13, units = "cm")


######
### ---- plot state-space dynamics for fig 1d----
###
x1 = 300
x2 = 580
y1 = .16
y2 = .3
y1bis = .05
xOffset = 25
yOffset= 15
colSquare = "yellow"
colLine = "white"
linSquare="dashed"
linLine = "dotted"
fontSize = 3
xLim = c(0,800)
yLim = c(0,.5)


npi.upper.plot = detectedPlot + scale_x_continuous(breaks = c(0, 400, 800), 
                                  labels = c("","",""),
                                  expand = c(0,0),
                                  limits = xLim) +
  scale_y_continuous(expand=c(0,0),limits=yLim)+
  ggtitle("") +
  labs(x="Strength of NPIs",y="Propotion infected")+
  
  theme(#axis.title.y=element_blank(),
        legend.key.size = unit(.5,"cm"),
        legend.text = element_text(size = 9),
        legend.position = "right")+
  annotate(geom="point",x=x1,y=y1,col=colSquare)+
  annotate(geom="point",x=x2,y=y1,col=colSquare)+
  annotate(geom="point",x=x1,y=y2,col=colSquare)+
  annotate(geom="point",x=x2,y=y2,col=colSquare) + 
  annotate(geom = "text",x=x2+xOffset,y1,label="italic(i)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x2+xOffset,y2,label="italic(iv)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x1-xOffset,y1,label="italic(ii)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x1-xOffset,y2,label="italic(iii)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom="segment",x=x2,xend = x1+25,y=y1,yend=y1,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type="closed"),linetype=linSquare) +
  annotate(geom="segment",x=x1,xend = x2-25,y=y2,yend=y2,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type = "closed"),linetype=linSquare)+
  annotate(geom="segment",x=x1,xend = x1,y=y1,yend=y2-.02,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type = "closed"),linetype=linSquare)+
  annotate(geom="point",x=x1,y=y1bis,col=colLine)+
  annotate(geom = "text",x=x1-xOffset,y1bis,label="italic(v)",size=fontSize,parse=T,hjust="outward",col=colLine)+
  annotate(geom="segment",x=x2,xend = x1+xOffset,y=y1,yend=y1bis,col=colLine,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type="closed"),linetype=linLine)


npi.upper.plot
ggsave("./plots/fig4_2021.png", npi.upper.plot, width = 6, height = 4, units = "in")



#npi.lower.plot = npi.upper.plot

#grid.arrange(npi.upper.plot,npi.lower.plot,ncol=1,left="Proportion infected",bottom="Strength of NPIs")


########################
#----  plot 1d phase space ----
########################

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
I50_w <- 400
maxDetected_logistic <- 0
Infected <- c(1:2000)

logisticR <- calculate_Rt(maxDetected=maxDetected_logistic, I50=I50,
                      maxLinks=maxLinks, Infected=Infected,
                      pCall=pCall, NDailyCalls=NDailyCalls,
                      pInfection=pInfection, Npop=Npop)
weakAlleeR <- calculate_Rt(maxDetected=maxDetected, I50=I50_w,
                      maxLinks=maxLinks, Infected=Infected,
                      pCall=pCall, NDailyCalls=NDailyCalls,
                      pInfection=pInfection, Npop=Npop)
alleeR <- calculate_Rt(maxDetected=maxDetected, I50=I50,
                      maxLinks=maxLinks, Infected=Infected,
                      pCall=pCall, NDailyCalls=NDailyCalls,
                      pInfection=pInfection, Npop=Npop)

infRdf <- data.frame(Infected = c(1:2000), logisticR = logisticR, alleeR = alleeR, weakAlleeR = weakAlleeR) %>%
  tidyr::pivot_longer(., cols = c("alleeR", "logisticR", "weakAlleeR"), names_to = "model",
                      values_to = "R") %>%
  dplyr::mutate(., model=factor(model, levels=c("logisticR", "weakAlleeR", "alleeR")))

# Find equilibrum points and put into dataframe to use in plot 
equilibriumPoints <- which(alleeR < 1)

point1 <- which(diff(equilibriumPoints) == max(diff(equilibriumPoints)))
point2 <- equilibriumPoints[point1+1]

equilibriumDf <- data.frame(Infected = c(point1, point2), R = c(1, 1),
                               type = c("unstable", "stable"))
arrowDf <- data.frame(Ii = c(point1-45, point1+120, Npop-100),
                      If = c(-60, point2-240, point2+100),
                      direction = c("left", "right", "left"),
                      R = 1.1)

thresholdText <- c("Epidemic \nthreshold")
immunityText <- c("Population \nimmunity")
#threshold_I_Text <- c("\nI*")

allee1D <- ggplot(infRdf, aes(x=Infected, y=R, color=model)) +
  geom_line(size=1) +
  scale_color_manual(values=c("#14b74b", "black", "#b33018"), name="",
                     labels=c("no NPI-Allee effect",
                                "weak NPI-Allee effect", "strong NPI-Allee effect")) +
  geom_hline(yintercept=1, size=1, linetype="dashed") +
  geom_segment(aes(x=point1, xend=point1, y=0, yend=2), color="red",
               linetype="dashed") +
  xlab("Proportion infected") +
  ylab(bquote('Reproductive number ' ~R[e] )) +
  geom_point(data=equilibriumDf, size=4, color="black") +
  geom_segment(data=dplyr::filter(arrowDf, direction == "right"),
               aes(x=Ii, xend=If, y=R, yend=R), color="#fa3d1b",
               arrow=arrow(length=unit(0.2, "cm")), size=0.9) +
  geom_segment(data=dplyr::filter(arrowDf, direction == "left"),
               aes(x=Ii, xend=If, y=R, yend=R), color= "#243faf",
               arrow=arrow(length=unit(0.2, "cm"), ends="last"), size=0.9) +
  annotate("text", x=point1+50, y=0.7, label=thresholdText, size=3.2, hjust="inward") +
  #annotate("text", x=point1+370, y=0.7, label=threshold_I_Text, size=3.2, hjust="inward",fontface="italic") +
  annotate("text", x=point2, y=0.7, label=immunityText, size=3.2, hjust="inward") +
  #annotate("text",x=point1+100,y=2.5,label="Allee effect \n (saturation of NPIs)",size=3,hjust="center") +
#  geom_text(aes(x = point2, y = 0.7), label = immunityText, color = "black",
 #           hjust = "inward", family = "sans", fontface = "plain" ) +
  theme_classic() +
  scale_x_continuous(breaks=c(0, 1000, 2000), labels=c(0, 0.5, 1)) + 
  theme(legend.position=c(.7,.8),text=element_text(size=12)) 

ggsave("./plots/allee1D.png", allee1D, width = 15, height = 12, units = "cm")


fig1_top_row = ggarrange(allee1D+theme(plot.margin=margin(8,8,8,8)),
                         space2Dplots, nrow=1,labels="AUTO")
ggsave("./plots/fig1.png", fig1_top_row, width=10, height=4, units="in")

