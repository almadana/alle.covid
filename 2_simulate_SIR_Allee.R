# Simulates a stochastic SIR model with and without Allee effect.
# Generates the plot from Figure 1C

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(segmented)
library(gridExtra)

# simulation parameters

set.seed(2691)
nRep <- 1000
minI <- 10
maxI <- 2000
popSizeStatsCounties <- c(mean = 5.2, sd = 0.5)
#popSizeStatsCountries <- c(mean = 7.1, sd = 0.5)
# epidemic values
p_se <- 0.8 # 1 es sin subexponencial
I50 <- 20
betaMax <- 0.5
gammaMax <- 4
# spread values
beta_dispersion <- 0.2
muImported <- 1
#p_se=0.8 # subexponencial
# sample populations
popSizesCounties <- round(10^rnorm(nRep, popSizeStatsCounties[1], popSizeStatsCounties[2]))
#popSizesCountries <- round(10^rnorm(nRep, popSizeStatsCountries[1], popSizeStatsCountries[2]))
minEpidemicLength <- 10

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
      I[t-1] <- 0
      I[t] <- 0
    }
    if(S[t-1] <= 0) {
      S[t-1] <- 0 
      S[t] <- 0
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
SIR_generator <- function(p_se, allee, popSizes, I50, betaMax = 1.4,
                          gammaMax = 5, dispersion = 0.2, muImported = 1) {
  outputDf <- data.frame()
  nRep <- length(popSizes)
  for (al in allee) {
    rep <- 1
    while (rep <= nRep) {
      N <- popSizes[rep]
      sim <- SIR_Allee(I0=10, betaMax=betaMax, gammaMax=gammaMax, p=p_se,
                                   durSim=200,  migration=muImported,
                                   dispersion=dispersion,
                                   I50=I50, Nsus=N, Allee = al)
      sim <- dplyr::mutate(sim, cumI = cumsum(newI + Imported),
                           cumImported = cumsum(Imported),
                           rep = rep, allee = al, p = p_se,
                           Population = N) %>%
        dplyr::filter(., cumI > minI & cumI < maxI)
      # crop up to day with largest daily count
      dailyCases <- diff(sim$cumI)
      maxDay <- which(dailyCases == max(dailyCases))[1]
      sim <- sim[1:maxDay,]
      if (nrow(sim) > 0) {
        sim$t <- c(1:nrow(sim))
        if (nrow(sim) > minEpidemicLength) {
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
    if.false <- F
    #sometime segmented returns an error, if case, run again
    tries <- 1
    while(if.false == F){
      tryCatch({
        fit <- segmented(lm(cumI ~ t), seg.Z = ~t, npsi = 1)
        print(tries)
        if.false <- T
      }, error = function(e){print(tries)
      }, finally = {print(tries)})
      tries <- tries + 1
      # if it can't fit the simulation, set everything to NA
      if (tries > 100) {
        fit <- lm(cumI ~ t)
      }
    }
    cfs <- summary(fit)
    slopeVals <- cfs$coefficients[,1]
    breakingPoint <- ifelse(is.null(cfs$psi), NA, cfs$psi[2])
    if (length(slopeVals) <= 2) {
      slopeVals <- c(slopeVals[1:2], NA)
    } else {
      # put final slope instead of difference in slopes
      slopeVals <- slopeVals[1:3]
      slopeVals[3] <- slopeVals[2] + slopeVals[3]
    }
    ff0 <- lm(cumI~ t)
    aic <- AIC(ff0,fit)[,2]
    wi <- exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
    # check if breaking point is significant by fitting lm
    # put together
    coefficients <- c(slopeVals, breakingPoint, wi[2])
    coefficients <- as.list(coefficients)
    names(coefficients) <- c("intercept", "slopeI", "slopeF",
                             "time.threshold","weighted.evidence")
  return(coefficients)
}

#### Run simulations and fit segmented for countie sims######
allee <- c(TRUE, FALSE)
simulationsCounties <- SIR_generator(p_se = p_se,
                             betaMax = betaMax,
                             gammaMax = gammaMax,
                             dispersion = beta_dispersion,
                             muImported = muImported,
                             allee = allee,
                             popSizes = popSizesCounties,
                             I50 = I50)

fitCoefsCounties <- group_by(simulationsCounties, allee, rep, p, Population) %>% 
  dplyr::summarise(., coefs = list(fit_segmented(cumI))) %>%
  ungroup(.) %>% 
  tidyr::unnest_wider(., coefs) %>% 
  dplyr::mutate(., slopeRatio = slopeF/slopeI, resta = slopeF - slopeI,
                angle = atan(abs((slopeF - slopeI)/(1-slopeF*slopeI))),
                Ithreshold = intercept + slopeI*time.threshold)

saveRDS(simulationsCounties, "./generated_data/SIR_dynamics_simulation.RDS") 
saveRDS(fitCoefsCounties, "./generated_data/SIR_dynamics_fit.RDS")

