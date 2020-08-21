library(segmented)
library(dplyr)
library(ggplot2)

# Countries Covid data
countriesCovid <- read.delim("./raw_data/casos.world.264.tab") %>%
  as_tibble(.)

minimumCumCases <- 100

#### Fit the segmented lines to the country data, in log(infected) vs log(time)
#### and analyze parameters and weight of evidence for breakpoint
#### Infected = Io·(time)^mu following Maiers et al. 2020
fit_country_segmented <- function(countriesCovid){
  #library("segmented")
  nCols <- dim(countriesCovid)[2]
  out <- NULL
  #par(mfrow=c(5,5), mar=c(2,2,2,0.5 ))
  for (country in 1:nrow(countriesCovid)){
    cumCases <- as.numeric(countriesCovid[country, 18:nCols]) #al 10 de julio
    # crop to beggining of epidemic
    epiInd <- which(cumCases >= 1)
    cumCases <- cumCases[epiInd]
    # get daily cases and day with max daily cases
    dailyCases <- diff(cumCases)
    maxDay <- which(dailyCases == max(dailyCases))[1]
    # select countries with more than 15 days of covid dynamics
    if (length(cumCases)>=15) {
      day10 <- which(cumCases > 10)[1]
      # Cut the initial 14 Days (many countries have weeks with few cases)
      initDay <- day10 + 14
      dynamicsLength <- maxDay - initDay + 1
      # original filter: if (length(14:id) > 15 & min(cumCases[14:id])>10 & max(cumCases)>200) {
      if (!is.na(day10) & dynamicsLength >= 15 & max(cumCases)>minimumCumCases) {
        croppedDyn <- cumCases[initDay:maxDay]
        logTime <- I(log10(1:length(croppedDyn)))
        logInfect <- log10(croppedDyn)
        fitSegmented <- segmented(obj = lm(logInfect ~ logTime), seg.Z = ~logTime, psi = 1)
        fit <- summary(fitSegmented)
        paramsSeg <- coefficients(fitSegmented)
        fitLinear <- lm(logInfect ~ logTime)
        aic <- AIC(fitLinear,fitSegmented)[,2]
        wi <- exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
        p.li <- coefficients(fitLinear)
        N.country <- countriesCovid[country,"Population"] # Population of the country
        # Estimate the number of active cases at the break point
        beforeCut <- which(logTime<fit$psi[2]) 
        if(length(beforeCut)>9) {
          # DH: ¿que hace esto?
          # If the time series before the break is larger than 10 days retain
          # the total infected 10 days before...
          # NOTA: Acá se estaba tomando la resta de los logaritmos como el
          # logaritmo de la resta me parece
          # I.active <- logInfect[max(beforeCut)] - logInfect[max(beforeCut)-9] # linea vieja
          I.active <- croppedDyn[max(beforeCut)] - croppedDyn[max(beforeCut)-9] # linea vieja
        } else if (length(beforeCut)==0) {
          I.active <- NA
        } else if (length(beforeCut)<=9) {
          # If the time series is shoreter take the total number of infected before the break
          I.active <- croppedDyn[max(beforeCut)]                        
        }
        if (length(fit$psi[2]) == 0) {
          timeThreshold <- NA
        } else {
          timeThreshold <- fit$psi[2]
        }
        countryId <- countriesCovid[country, "UID"]
        countryName <- countriesCovid[country, "Country_Region"]
        countryVec <- data.frame(countryId, countryName, N.country[[1]], max(cumCases),
                       paramsSeg[1], paramsSeg[2], paramsSeg[3],
                       timeThreshold, wi[2],
                       fit$r.squared, I.active)
        colnames(countryVec) <- c("Id", "Country", "Population", "I.max",
                                  "intercept", "slopeI", "slopeF", "time.threshold",
                                  "weighted.evidence","R.sqrt", "I.active")
        out <- rbind(out, countryVec)
        fitLinear <- NULL
        fitSegmented <- NULL
      }
    }      
  }
  # calculate infected number at the time of breakout
  intercepts <- as.numeric(out[,"intercept"])
  slope <- as.numeric(out[,"slopeI"])
  breakTime <- as.numeric(out[,"time.threshold"])
  Inf.total <- 10^(intercepts + slope * breakTime)
  out <- cbind(out, Inf.total)
  # the segmented package returns difference in slopes. Turn to final slope
  out$slopeF <- out$slopeF + out$slopeI
  # name the columns
  rownames(out) <- NULL
  return(out)
}      

### Fit the countries and save df
countrySegmentedFit <- fit_country_segmented(countriesCovid) %>%
  dplyr::filter(., !is.na(weighted.evidence))

saveRDS(countrySegmentedFit, "./generated_data/countriesFit.RDS")

