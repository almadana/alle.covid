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
library(ggpubr)
source("./functions_static_model.R")


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
I50 <- 200

# ranges for 
maxDetected_Vec <- seq(0, 800, 5)
dailyCalls_Vec <- seq(20, 1500, 5)
maxLinks_Vec <- seq(2, 30, 0.1)
pInfection_Vec <- seq(0.1, 0.3, 0.0015)

# Space plot for maximum number of people that can be detected
detectedLong <- Rt_detected(maxLinks=maxLinks, maxDetected=maxDetected_Vec,
                            I50=I50, pCall=pCall,
                            pInfection=pInfection, NDailyCalls=NDailyCalls,
                            Npop=Npop, propInfectedMax=propInfectedMax) %>%
                enlongate_matrix(., c("Infected", "maxDetected"), infectedSubsample)

detectedPlot <- dplyr::rename(detectedLong, measure=maxDetected) %>%
  dplyr::mutate(Infected=Infected/Npop) %>% 
  plot_phase_space(.) +
  ggtitle("Tracing capacity") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.title.x = element_text(size = 10)) +
  scale_x_continuous(breaks = c(0, 400, 800),
                     name = bquote("Detection capacity (K)"),
                     expand = c(0,0))

# Space plot for maximum number of calls (speed of detection)
callsLong <- Rt_calls(maxLinks=maxLinks, maxDetected=maxDetected,
                      I50=I50, pCall=pCall, pInfection=pInfection,
                      NDailyCalls=dailyCalls_Vec, Npop=Npop,
                      propInfectedMax=propInfectedMax) %>%
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
linksLong <- Rt_links(maxLinks=maxLinks_Vec, maxDetected=maxDetected,
                      I50=I50, pCall=pCall, pInfection=pInfection,
                      NDailyCalls=NDailyCalls, Npop=Npop,
                      propInfectedMax=propInfectedMax) %>% 
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
infectionLong <- Rt_pInfection(maxLinks=maxLinks, maxDetected=maxDetected,
                               I50=I50, pCall=pCall,
                               pInfection=pInfection_Vec,
                               NDailyCalls=NDailyCalls, Npop=Npop,
                               propInfectedMax=propInfectedMax) %>%
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


ggsave("./plots/phase_space.pdf", space2Dplots, width = 15, height = 13, units = "cm")


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
ggsave("./plots/fig4_2021.pdf", npi.upper.plot, width = 6, height = 4, units = "in")



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
# Df for 3 arrows
#arrowDf <- data.frame(Ii = c(point1-45, point1+120, Npop-100),
#                      If = c(-60, point2-240, point2+100),
#                      direction = c("left", "right", "left"),
#                      R = 1.1)
# Df for 2 arrows
arrowDf <- data.frame(Ii = c(point1-45, point1+45),
                      If = c(point1-150, point1+150),
                      direction = c("left", "right"),
                      R = 1.07)


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

ggsave("./plots/allee1D.pdf", allee1D, width = 15, height = 12, units = "cm")


fig1_top_row = ggpubr::ggarrange(allee1D+theme(plot.margin=margin(8,8,8,8)),
                         space2Dplots, nrow=1,labels="AUTO")
ggsave("./plots/fig1.pdf", fig1_top_row, width=10, height=4, units="in")



#####################
# Get the proportion of undetected infections and
# the time to testing changes as functions of I
#####################

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

alleeR <- calculate_Rt(maxDetected=maxDetected, I50=I50,
                      maxLinks=maxLinks, Infected=Infected,
                      pCall=pCall, NDailyCalls=NDailyCalls,
                      pInfection=pInfection, Npop=Npop)

propFound <- calculate_proportion_found(maxDetected=maxDetected,
                                        I50=I50, Infected=Infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls)

meanDay <- calculate_detection_time(maxDetected=maxDetected,
                                        I50=I50, Infected=Infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls)

foundDf <- data.frame(Infected=Infected, R=alleeR, 
                      propFound=propFound, notFound=Infected-Infected*propFound,
                      time=meanDay)

equilibriumPoints <- which(alleeR < 1)
point1 <- which(diff(equilibriumPoints) == max(diff(equilibriumPoints)))
point2 <- equilibriumPoints[point1+1]

foundPlot <- ggplot(foundDf, aes(x=Infected, y=1-propFound)) +
  geom_line(aes(color="Undetected infections")) +
  geom_line(aes(y=R/max(alleeR), color="Re")) +
  geom_hline(yintercept=1/max(alleeR), linetype="dashed", color="red") +
  geom_vline(xintercept=point1, linetype="dashed", color="red") +
#  ylab() +
  scale_y_continuous(name="Proportion of undetected infections",
                     sec.axis=sec_axis(~.*max(alleeR), name="Re")) +
  scale_color_manual(name=element_blank(), values=c("red", "black")) +
  theme_bw() +
  theme(axis.title.y.right=element_text(color="red"),
        axis.text.y.right=element_text(color="red"),
        axis.line.y.right=element_line(colour="red"),
        axis.ticks.y.right=element_line(colour="red"),
        legend.position="top")

timeFoundPlot <- ggplot(foundDf, aes(x=Infected, y=time)) +
  geom_line(aes(color="Time of contact")) +
  geom_line(aes(y=R/max(alleeR)*max(meanDay), color="Re")) +
  geom_hline(yintercept=1/max(alleeR)*max(meanDay), linetype="dashed", color="red") +
  geom_vline(xintercept=point1, linetype="dashed", color="red") +
#  ylab() +
  scale_y_continuous(name="Mean time to successful contact (days)",
                     sec.axis=sec_axis(~.*max(alleeR)/max(meanDay), name="Re")) +
  scale_color_manual(name=element_blank(), values=c("red", "black")) +
  theme_bw() +
  theme(axis.title.y.right=element_text(color="red"),
        axis.text.y.right=element_text(color="red"),
        axis.line.y.right=element_line(colour="red"),
        axis.ticks.y.right=element_line(colour="red"),
        legend.position="top")

indicatorsPlot <- ggpubr::ggarrange(foundPlot, timeFoundPlot, nrow=1, ncol=2,
labels=c("A", "B"))

ggsave("./plots/supp2_indicadores_outbreak_estatico.pdf", indicatorsPlot,
       width=18, height=8, units="cm")

#############################
# Same for weak allee
#############################
weakAlleeR <- calculate_Rt(maxDetected=maxDetected, I50=I50_w,
                      maxLinks=maxLinks, Infected=Infected,
                      pCall=pCall, NDailyCalls=NDailyCalls,
                      pInfection=pInfection, Npop=Npop)

propFound <- calculate_proportion_found(maxDetected=maxDetected,
                                        I50=I50_w, Infected=Infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls)

meanDay <- calculate_detection_time(maxDetected=maxDetected,
                                        I50=I50_w, Infected=Infected,
                                        pCall=pCall,
                                        NDailyCalls=NDailyCalls)

foundDf <- data.frame(Infected=Infected, R=weakAlleeR, 
                      propFound=propFound, notFound=Infected-Infected*propFound,
                      time=meanDay)


foundPlot <- ggplot(foundDf, aes(x=Infected, y=1-propFound)) +
  geom_line(aes(color="Undetected infections")) +
  geom_line(aes(y=R/max(alleeR), color="Re")) +
#  geom_hline(yintercept=1/max(alleeR), linetype="dashed", color="red") +
#  geom_vline(xintercept=point1, linetype="dashed", color="red") +
#  ylab() +
  scale_y_continuous(name="Proportion of undetected infections",
                     sec.axis=sec_axis(~.*max(alleeR), name="Re")) +
  scale_color_manual(name=element_blank(), values=c("red", "black")) +
  theme_bw() +
  theme(axis.title.y.right=element_text(color="red"),
        axis.text.y.right=element_text(color="red"),
        axis.line.y.right=element_line(colour="red"),
        axis.ticks.y.right=element_line(colour="red"),
        legend.position="top")

timeFoundPlot <- ggplot(foundDf, aes(x=Infected, y=time)) +
  geom_line(aes(color="Time of contact")) +
  geom_line(aes(y=R/max(alleeR)*max(meanDay), color="Re")) +
#  geom_hline(yintercept=1/max(alleeR)*max(meanDay), linetype="dashed", color="red") +
#  geom_vline(xintercept=point1, linetype="dashed", color="red") +
#  ylab() +
  scale_y_continuous(name="Mean time to successful contact (days)",
                     sec.axis=sec_axis(~.*max(alleeR)/max(meanDay), name="Re")) +
  scale_color_manual(name=element_blank(), values=c("red", "black")) +
  theme_bw() +
  theme(axis.title.y.right=element_text(color="red"),
        axis.text.y.right=element_text(color="red"),
        axis.line.y.right=element_line(colour="red"),
        axis.ticks.y.right=element_line(colour="red"),
        legend.position="top")

indicatorsPlot <- ggpubr::ggarrange(foundPlot, timeFoundPlot, nrow=1, ncol=2)

ggsave("./plots/supp2_indicadores_outbreak_estatico_weakAllee.pdf", indicatorsPlot,
       width=18, height=8, units="cm")

