library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
source("./functions_static_model.R")


###########################
# Reproduce figure 1: R0 vs Infected individuals
###########################

Npop <- 2000
I0_vec <- c(1:2000)

# Logistic model, no variation in beta/gamma
beta <- 0.7
gamma <- 4
S0_vec <- Npop-I0_vec
logisticR <- beta*gamma*(S0_vec/Npop)

# Strong allee effect
betaMax <- 0.7
I50b <- 100
gammaMax <- 4
I50g <- 100
alleeR <- NULL
for (ind in c(1:length(I0_vec))) {
  I0 <- I0_vec[ind]
  S0 <- S0_vec[ind]
  beta <- betaMax*I0/(I50b+I0)
  gamma <- gammaMax*I0/(I50g+I0)
  alleeR[ind] <- beta*gamma*(S0/Npop)
}


# Weak allee effect
betaMax <- 0.7
I50b <- 100
gammaMax <- 4
I50g <- 100
weakAlleeR <- NULL
betaMin <- 0.4
gammaMin <- 3
for (ind in c(1:length(I0_vec))) {
  I0 <- I0_vec[ind]
  S0 <- S0_vec[ind]
  beta <- betaMin + (betaMax-betaMin)*I0/(I50b+I0)
  gamma <- gammaMin + (gammaMax-gammaMin)*I0/(I50g+I0)
  weakAlleeR[ind] <- beta*gamma*(S0/Npop)
}


infRdf <- data.frame(Infected=I0_vec, logisticR=logisticR, alleeR=alleeR,
                     weakAlleeR=weakAlleeR) %>%
  tidyr::pivot_longer(., cols = c("alleeR", "logisticR", "weakAlleeR"), names_to = "model",
                      values_to = "R") %>%
  dplyr::mutate(., model=factor(model, levels=c("logisticR", "weakAlleeR", "alleeR")))


# Find equilibrum points and put into dataframe to use in plot 
equilibriumPoints <- which(alleeR < 1)

point1 <- which(diff(equilibriumPoints) == max(diff(equilibriumPoints)))
point2 <- equilibriumPoints[point1+1]

equilibriumDf <- data.frame(Infected = c(point1, point2), R = c(1, 1),
                               type = c("unstable", "stable"))
# 2 arrow df
#arrowDf <- data.frame(Ii = c(point1-45, point1+45, Npop-100),
#                      If = c(-60, point2-240, point2+100),
#                      direction = c("left", "right", "left"),
#                      R = 1.1)
# 3 arrow df
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
                     labels=c("SIR without NPIs",
                                "SIR with saturable NPIs, weak Allee",
                                "SIR with saturable NPIs, strong Allee")) +
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



###########################
# Reproduce figure 2: R0 vs Infected individuals
###########################

# general parameters
betaMax <- 0.7
I50b <- 100
gammaMax <- 4
I50g <- 100
alleeR <- NULL

# betaMax sweep
betaMax_vec <- seq(0.3, 1, 0.002)
betaMaxDf <- data.frame()
for (bm in c(1:length(betaMax_vec))){ 
  betaMax_sample <- betaMax_vec[bm]
  for (ind in c(1:length(I0_vec))){
    I0 <- I0_vec[ind]
    S0 <- S0_vec[ind]
    beta <- betaMax_sample*I0/(I50b+I0)
    gamma <- gammaMax*I0/(I50g+I0)
    alleeR[ind] <- beta*gamma*(S0/Npop)
  }
  tempDf <- data.frame(measure=betaMax_sample, Infected=I0_vec, Rt_exp=alleeR)
  betaMaxDf <- rbind(betaMaxDf, tempDf)
}
betaMaxDf <- dplyr::mutate(betaMaxDf, Rt_exp=Rt_exp>1, Infected=Infected/Npop)
betaMaxPlot <- plot_phase_space(inputDF=betaMaxDf) +
  theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
        axis.title.x=element_text(size=10)) +
  scale_x_continuous(name=expression(beta[max]),
                     expand=c(0,0))

# gammaMax sweep
gammaMax_vec <- seq(2, 6, 0.01)
gammaMaxDf <- data.frame()
for (bm in c(1:length(gammaMax_vec))){ 
  gammaMax_sample <- gammaMax_vec[bm]
  for (ind in c(1:length(I0_vec))){
    I0 <- I0_vec[ind]
    S0 <- S0_vec[ind]
    beta <- betaMax*I0/(I50b+I0)
    gamma <- gammaMax_sample*I0/(I50g+I0)
    alleeR[ind] <- beta*gamma*(S0/Npop)
  }
  tempDf <- data.frame(measure=gammaMax_sample, Infected=I0_vec, Rt_exp=alleeR)
  gammaMaxDf <- rbind(gammaMaxDf, tempDf)
}
gammaMaxDf <- dplyr::mutate(gammaMaxDf, Rt_exp=Rt_exp>1, Infected=Infected/Npop)
gammaMaxPlot <- plot_phase_space(inputDF=gammaMaxDf) +
  theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
        axis.title.x=element_text(size=10)) +
  scale_x_continuous(name=expression(tau[max]),
                     expand=c(0,0))

# I50b sweep
I50b_vec <- seq(20, 500, 2)
betaI50Df <- data.frame()
for (bm in c(1:length(I50b_vec))){ 
  I50b_sample <- I50b_vec[bm]
  for (ind in c(1:length(I0_vec))){
    I0 <- I0_vec[ind]
    S0 <- S0_vec[ind]
    beta <- betaMax*I0/(I50b_sample+I0)
    gamma <- gammaMax*I0/(I50g+I0)
    alleeR[ind] <- beta*gamma*(S0/Npop)
  }
  tempDf <- data.frame(measure=I50b_sample, Infected=I0_vec, Rt_exp=alleeR)
  betaI50Df <- rbind(betaI50Df, tempDf)
}
betaI50Df <- dplyr::mutate(betaI50Df, Rt_exp=Rt_exp>1, Infected=Infected/Npop)
betaI50Plot <- plot_phase_space(inputDF=betaI50Df) +
  theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
        axis.title.x=element_text(size=10)) +
  scale_x_continuous(name=expression(I[50][beta]),
                     expand=c(0,0))


# I50b sweep
I50g_vec <- seq(20, 500, 4)
gammaI50Df <- data.frame()
for (bm in c(1:length(I50g_vec))){ 
  I50g_sample <- I50g_vec[bm]
  for (ind in c(1:length(I0_vec))){
    I0 <- I0_vec[ind]
    S0 <- S0_vec[ind]
    beta <- betaMax*I0/(I50b+I0)
    gamma <- gammaMax*I0/(I50g_sample+I0)
    alleeR[ind] <- beta*gamma*(S0/Npop)
  }
  tempDf <- data.frame(measure=I50g_sample, Infected=I0_vec, Rt_exp=alleeR)
  gammaI50Df <- rbind(gammaI50Df, tempDf)
}
gammaI50Df <- dplyr::mutate(gammaI50Df, Rt_exp=Rt_exp>1, Infected=Infected/Npop)
gammaI50Plot <- plot_phase_space(inputDF=gammaI50Df) +
  theme(plot.title=element_text(hjust=0.5, size=11, face="bold"),
        axis.title.x=element_text(size=10)) +
  scale_x_continuous(name=expression(I[50][tau]),
                     expand=c(0,0))

space2Dplots <- grid.arrange(
                     arrangeGrob(betaMaxPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 gammaMaxPlot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                  betaI50Plot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 gammaI50Plot + theme(legend.position="none",
                                                      axis.title.y=element_blank(),
                                                      plot.margin=margin(8,8,8,8)),
                                 ncol=2, left = "Proportion infected")#,legend,
                     #heights=c(10, .5)
                     )



fig1_reproduced_SIR <- ggpubr::ggarrange(allee1D+theme(plot.margin=margin(8,8,8,8)),
                         space2Dplots, nrow=1,labels="AUTO")
ggsave("./plots/fig1_SIR_version.pdf", fig1_reproduced_SIR, width=10, height=4, units="in")

ggsave("./plots/phase_space_SIR.pdf", space2Dplots, width=10, height=7, units="in")
ggsave("./plots/allee1D_SIR.pdf", allee1D, width=10, height=7, units="in")


