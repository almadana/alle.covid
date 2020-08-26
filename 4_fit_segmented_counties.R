#### Dara bases from The New York Times and from US Census
library(segmented)
library(dplyr)
library(ggplot2)

set.seed(2691)
minimumCases <- 200

# US counties COVID data
countiesCovid <- read.csv("./raw_data/us-counties.csv", colClasses = c("fips" = "character")) %>%
  as_tibble(.) %>%
  dplyr::filter(., county != "Unknown")

# Demographic data and counties identification
countiesDem <- read.csv("./raw_data/co-est2019-alldata.csv", header=TRUE,
                        colClasses = c("STATE" = "character", "COUNTY" = "character")) %>%
  as_tibble(.)

# make counties Id for the demographics df
countiesId <- paste0(countiesDem$STATE, countiesDem$COUNTY)
countyPop <- as_tibble(cbind(countiesId, countiesDem[,6:8])) # id of counties
colnames(countyPop)[1] <- "fips"

#### Fit the segmented lines to the county data, in log(infected) vs log(time)
#### and analyze parameters and weight of evidence for breakpoint
#### Infected = Io·(time)^mu following Maiers et al. 2020
fit_county_segmented <- function(countyPop, countiesCovid){
  #library("segmented")
  out <- NULL
  #par(mfrow=c(5,5), mar=c(2,2,2,0.5 ))
  for (countyId in countyPop$fips){
    countyId <- as.character(countyId)
    dynInd <- which(countiesCovid$fips == countyId)
    # select counties with more than 15 days of covid dynamics
    if (length(dynInd)>=15) {
      cumCases <- dplyr::pull(countiesCovid[dynInd,5])
      dailyCases <- diff(cumCases)
      day10 <- which(cumCases > 10)[1]
      # Cut the initial 14 Days (many counties have weeks with few cases)
      initDay <- day10 + 14
      maxDay <- which(dailyCases == max(dailyCases))[1]
      dynamicsLength <- maxDay - initDay + 1
      # original filter: if (length(14:id) > 15 & min(cumCases[14:id])>10 & max(cumCases)>200) {
      if (!is.na(day10) & dynamicsLength >= 15 & max(cumCases)>minimumCases) {
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
        N.county <- countyPop[which(countyPop$fips == countyId),4] # Population of the county
        # Estimate the number of active cases at the break point
        beforeCut <- which(logTime<fit$psi[2]) 
        if(length(beforeCut)>9) {
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
        countyVec <- data.frame(countyId, N.county[[1]], max(cumCases),
                       paramsSeg[1], paramsSeg[2], paramsSeg[3],
                       timeThreshold, wi[2],
                       fit$r.squared, I.active)
        colnames(countyVec) <- c("Id", "Population", "I.max", "intercept",
                           "slopeI", "slopeF", "time.threshold",
                           "weighted.evidence","R.sqrt", "I.active")
        out <- rbind(out, countyVec)
        fitLinear <- NULL
        fitSegmented <- NULL
      }
    }      
  }
  # calculate infected number at the time of breakout
  intercepts <- as.numeric(out[,4])
  slope <- as.numeric(out[,5])
  breakTime <- as.numeric(out[,7])
  Inf.total <- 10^(intercepts + slope * breakTime)
  out <- cbind(out, Inf.total)
  # the segmented package returns difference in slopes. Turn to final slope
  out$slopeF <- out$slopeF + out$slopeI
  # name the columns
  rownames(out) <- NULL
  return(out)
}      

### Fit the counties and save df
countySegmentedFit <- fit_county_segmented(countyPop = countyPop,
                                           countiesCovid = countiesCovid) %>%
  dplyr::filter(., !is.na(weighted.evidence))

saveRDS(countySegmentedFit, "./generated_data/countiesFit.RDS")

