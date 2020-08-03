e##################################################################
##### Dynamic

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(segmented)
library(gridExtra)

# SIR simulation with Allee effect
SIR_Allee <- function(I0, betaMax=1.4, gammaMax=5, p, durSim,
                                  Allee = TRUE, migration = 1,
                                  dispersion = 0.2, I50 = 10, Nsus = 100000){
  S <- rep(NA, durSim)
  I <- rep(NA, durSim)
  R <- rep(NA, durSim)
  newI <- rep(NA, durSim)
  Imported <- rep(NA, durSim)
  # assign initial values
  S[1] <- Nsus
  I[1] <- I0
  R[1] <- 0
  newI[1] <- 0
  Imported[1] <- 0
  Npop <- Nsus + I0
  # preassign params for non-Allee
  beta <- betaMax
  gamma <- gammaMax
  for(t in 2:durSim){
    if (I[t-1] <= 0) {
      I[t-1] <- I[t] <- 0
    }
    if(S[t-1] <= 0) {
      S[t-1] <- S[t] <- 0
    }
    if (Allee == TRUE) {
      # adjust beta and gamma to infected levels if Allee
      beta <- betaMax * I[t-1] / (I[t-1] + I50)
      gamma <- gammaMax * I[t-1] / (I[t-1] + I50)
    }
    # get expected changes 
    newIExpected <- beta * S[t-1] * (I[t-1])^p / Npop
    recoveredIExpected <- (1/gamma) * I[t-1]
    # sample the actual changes from random variables
    newISample <- rpois(n = 1, lambda  = newIExpected)
    newI[t] <- min(newISample, S[t-1]) # ensure no more new infections than susceptible
    recoveredI <- rnbinom(n = 1, size = dispersion, mu = recoveredIExpected)
    recoveredI <- min(I[t-1], recoveredI, na.rm=T) # no more recovered than infected
    Imported[t] <- rnbinom(n = 1, size = dispersion, mu = migration)
    #new.I<-rnbinom(n = 1, size = dispersion, mu = new.I)
    S[t] <- S[t-1] - newI[t]
    I[t] <- I[t-1] + newI[t] + Imported[t] - recoveredI
    R[t] <- R[t-1] + recoveredI
    #out[i,]<-c(out[t-1,1]-new.I, out[t-1,2]+new.I-recoveredI+imported,out[t-1,3]+recoveredI, new.I, imported)
    }
  # store results in dataframe
  out <- data.frame(t = c(1:length(S)), S = S, I = I, R = R,
                    newI = newI, Imported = Imported)
  return(out)
}

# generate repetitions of simulations
SIR_generator <- function(p_se, allee, nRep, I50) {
  
  #N=round(10^rnorm(nRepeats,pob_paises[1],pob_paises[2]))  #valores para paises
  
  
  outputDf <- data.frame()
  for (al in allee) {
    rep <- 1
    while (rep <= nRep) {
      #N=round(10^rnorm(1,pob_condados[1],pob_condados[2]))  #valores para condados
      N=100000
      #nRows=length(alleeCoefs)
      #para explorar
      #subExpCoefs = runif(1)*0.2 + 0.7 #uniforme entre 0.7 y 0.9
      subExpCoefs=0.8
      #betamax = runif(1)*0.11 + 1.07
      betamax = 1.12
      #gammamax = runif(1)*3 + 4
      gammamax=5
      
      sim <- SIR_Allee(I0=10, betaMax=betamax, gammaMax=gammamax , p=subExpCoefs,
                                   durSim=200,  migration=1, dispersion=0.2,
                                   I50=I50, Nsus=N, Allee = al)
      sim <- dplyr::mutate(sim, cumI = cumsum(newI),
                           cumImported = cumsum(Imported),
                           rep = rep, allee = al, p = p_se) %>%
        dplyr::filter(., cumI > minI & cumI < maxI)
      if (nrow(sim) > 0) {
        sim$t <- c(1:nrow(sim))
        if (nrow(sim) > 10) {
          outputDf <- rbind(outputDf, sim)
          rep <- rep + 1
        }
      }
    }
  }
  return(outputDf)
}

# fit segmented lm to dynamics of single run
fit_segmented <- function(cumI) {
    # fit segmented model  
    cumI <- log10(cumI)
    #t <- log10(c(1:length(cumI)))
    t <- log10(c(1:length(cumI)))
    fit <- segmented(lm(cumI ~ t), seg.Z = ~t, npsi = 1)
    cfs <- summary(fit)
    slopeVals <- cfs$coefficients[,1]
    if (length(slopeVals) <= 2) {
      slopeVals <- c(slopeVals[1:2], NA)
    } else {
      slopeVals <- slopeVals[1:3]
      slopeVals[3] <- slopeVals[2] + slopeVals[3]
    }
    breakingPoint <- ifelse(is.null(cfs$psi), NA, cfs$psi[2])
    
    ff0<-lm(cumI~ t)
    aic<-AIC(ff0,fit)[,2]
    wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
    
    # check if breaking point is significant by fitting lm
    #lmFit <- lm(cumI ~ t)
    #pValue <- davies.test(lmFit)$p.value
    # put together
    coefficients <- c(slopeVals, breakingPoint,wi[2],fit)
    coefficients <- as.list(coefficients)
    names(coefficients) <- c("intercept", "slopeI", "slopeF",
                             "breakPoint","weighted.evidence")
  return(coefficients)
}

#fit_segmented <- function(cumI) {
#    # fit segmented model  
#    cumI <- log10(cumI)
#    #t <- log10(c(1:length(cumI)))
#    t <- c(1:length(cumI))
#    fit <- segmented(lm(cumI ~ t), seg.Z = ~t, npsi = 1)
#    cfs <- summary(fit)
#    slopeVals <- slope(fit)
#    slopeVals <- slopeVals$t[,1]
#    if (length(slopeVals) <= 1) {
#      slopeVals <- c(slopeVals, NA)
#    }
#    intercept <- coef(cfs)[1,1]
#    breakingPoint <- ifelse(is.null(cfs$psi), NA, cfs$psi[2])
#    coefficients <- c(intercept, slopeVals, breakingPoint)
#    # check if breaking point is significant by fitting lm
#    #lmFit <- lm(cumI ~ t)
#    #pValue <- davies.test(lmFit)$p.value
#    # put together
#    #coefficients <- as.list(c(coefficients, pValue))
#    coefficients <- as.list(coefficients)
#    names(coefficients) <- c("intercept", "slopeI", "slopeF",
#                             "breakPoint")
#  return(coefficients)
#}


# simulation parameters

set.seed(2691)
nRep <- 4000
minI <- 10
maxI <- 1000
p_se <- 0.8 # 1 es sin subexponencial
I50 <- 10
nPlotDyn <- 20
#p_se=0.8 # subexponencial


#### Run simulations and plot data ######
allee <- c(TRUE, FALSE)


pob_condados = c(4.49,0.68)
pob_paises = c(6.66,1.07)

simulations <- SIR_generator(p_se = p_se, allee = allee, nRep = nRep, I50 = I50)

#cumI <- filter(simulations, allee == TRUE & rep == 1)[["cumI"]]
#alSim <- filter(simulations, allee == TRUE)
#length(unique(alSim$rep))
#for (i in c(1:100)) {
#  cumI <- filter(simulations, allee == FALSE & rep == i)[["cumI"]]
#  k <- fit_segmented(cumI)
#}


fitCoefs <- 
group_by(simulations, allee, rep, p) %>% 
  dplyr::summarise(., coefs = list(fit_segmented(cumI))) %>%
  ungroup(.) %>% 
  tidyr::unnest_wider(., coefs) %>% 
  dplyr::mutate(., cociente = slopeF/slopeI, resta = slopeF - slopeI,
                angle = atan(abs((slopeF - slopeI)/(1-slopeF*slopeI))),
                Ithreshold = intercept + slopeI*breakPoint)
#  dplyr::filter(., cociente > 0)

fitCoefs$allee.f=as.factor(fitCoefs$allee)
levels(fitCoefs$allee.f)=c("without Allee effect","with Allee effect")
fitCoefs %>% group_by(allee) %>% summarize_at(vars(weighted.evidence),.funs = list(medianWE=median,meanWE=mean,minWE=min,sdWE=sd))


plotChange <- fitCoefs %>%
  dplyr::mutate(., allee = c("w/o Allee", "Allee")[as.integer(allee)+1]) %>%
  ggscatterhist(., x = "slopeI", y = "slopeF", color = "allee",
                alpha = 0.5, size = 3, margin.plot = "boxplot",
                palette = c("#b33018", "#14b74b"),
                ggtheme = theme_bw(), xlab="Initial slope",
                ylab="Slope after threshold",
                legend.title = "",
                margin.params = list(color = "allee"))


sampleRepetitions <- sample(unique(simulations$rep), nPlotDyn)

sample4 = sampleRepetitions[order(sampleRepetitions)]
sample4 = sample4[c(1,2,5,6)]

Nmax=100000

simulations_fit = merge(simulations,fitCoefs,by=c("rep","allee"))

simulations_fit = simulations_fit %>% mutate(t_bp = ceiling(10^breakPoint),cumI_fit = ifelse(t<10^breakPoint, 10^(intercept + log10(t)*slopeI),10^(intercept - log10(t_bp)*slopeF +log10(t_bp-1)*slopeI +log10(t)*slopeF ))) 

#plotsDynAllee <- dplyr::filter(simulations_fit, allee == TRUE & rep %in% sampleRepetitions) %>%
plotsDynAllee <- dplyr::filter(simulations_fit, allee == TRUE & rep %in%sample4) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#b33018") +
#  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  #ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("with Allee effect") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="pink")

plotsDynAllee
  
#plotsDynNonAllee <- dplyr::filter(simulations_fit, allee == FALSE & rep %in% sampleRepetitions) %>%
plotsDynNonAllee <- dplyr::filter(simulations_fit, allee == FALSE & rep %in% sample4) %>% View()
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#14b74b") +
#  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  xlab("time (days)") +
  #ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("without Allee effect") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="green")

dynsPlots <- grid.arrange(plotsDynAllee, plotsDynNonAllee, ncol=2,
                          left = "Cumulative infected")

ggsave("simDynamics.png", dynsPlots, width = 18, height = 11, units = "cm")
ggsave("simDynamics.pdf", dynsPlots, width = 18, height = 11, units = "cm")
#ggsave("slopeAnalysis.png", width = 12, height = 12, units = "cm")


