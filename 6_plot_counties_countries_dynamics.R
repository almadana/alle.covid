library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

fileAllDyns <- "./plots/data_dynamics.png"

fileSuppDyns <- "./plots/all_counties_dyn.png"

nExamples <- 8 # cuantos de cada
######################################
####### Plot counties dynamics #######
######################################

# load data
countiesFit <- readRDS("./generated_data/countiesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI) %>%
  as_tibble(.)

fittedCounties <- as.character(countiesFit$Id)

countiesCovid <- read.csv("./raw_data/us-counties.csv",
                          colClasses = c("fips" = "character")) %>%
  as_tibble(.) %>%
  dplyr::filter(., county != "Unknown" & cases > 10) %>%
  dplyr::group_by(., fips) %>% 
  dplyr::mutate(., t = row_number(), peak = which.max(diff(cases))[1]) %>% 
  dplyr::filter(., t > 14, t <= peak) %>% 
  dplyr::mutate(., t = t - 14) %>%
  dplyr::rename(., Id = fips) %>%
  dplyr::as_tibble(.)

stateAbrev <- read.csv("./raw_data/county_abbreviations.csv")

# add extended county name
countiesCovid <- merge(countiesCovid, stateAbrev, by = "state", all.x = T) %>%
  dplyr::mutate(., county.full = paste(county, code, sep=", ")) %>%
  dplyr::filter(., Id %in% fittedCounties)

countiesAll <- merge(countiesCovid, countiesFit, by="Id")

# calculate the fitted values for the countie dynamics
countiesAll <- countiesAll %>%
  dplyr::mutate(., t_bp = ceiling(10^time.threshold),
   cumI_fit = ifelse(t<10^time.threshold,
                     10^(intercept + log10(t)*slopeI),
                     10^(intercept - log10(t_bp)*slopeF +
                         log10(t_bp) * slopeI + log10(t) * slopeF))) 

# get counties with the least and most pronounced breakpoints
increasingRatio <- dplyr::filter(countiesFit, slopeRatio > 1 &
                                 Id != "13087")
slopeOrder <- base::order(increasingRatio$slopeRatio)
nCounties <- length(slopeOrder)
#ind <- round(c(seq(1, nCounties, nCounties/(nExamples-1)), nCounties))
ind <- c(seq(1, nExamples/2, 1), seq(nCounties-nExamples/2+1, nCounties, 1))
examplesIndices <- slopeOrder[ind]
examplesId <- as.character(increasingRatio$Id[examplesIndices])

plotsCounties <- dplyr::filter(countiesAll, Id %in% examplesId) %>%
  ggplot(., aes(x = t, y = cases)) +
  geom_point(color = "#000080") +
  facet_wrap(~county.full, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("time (days)") +
  #ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("U.S. counties") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

plotsCounties

######################################
####### Plot countries dynamics #######
######################################

countriesFit <- readRDS("./generated_data/countriesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI, Country = as.character(Country)) %>%
  dplyr::mutate(., Country = ifelse(Country %in% "Korea, South","South Korea",Country)) %>%
  as_tibble(.)

fittedCountries <- as.character(countriesFit$Id)

countriesCovid <- read.delim("./raw_data/casos.world.264.tab") %>%
  dplyr::rename(., Id = UID) %>%
  dplyr::filter(., Id %in% fittedCountries) %>%
  dplyr::as_tibble(.) %>%
  dplyr::select(., -X) %>%
  dplyr::select(., Id, Province_State, Country_Region, starts_with("X")) %>% 
  tidyr::pivot_longer(., starts_with("X")) %>%
  group_by(., Province_State, Country_Region) %>%
  dplyr::rename(., cumI = value) %>%
  dplyr::filter(., cumI > 10) %>%
  dplyr::mutate(., t = 1:n(), peak = which.max(diff(cumI))[1]) %>%
  ungroup(.) %>%
  dplyr::filter(., t > 14, t <= peak) %>% 
  dplyr::mutate(., t = t - 14) %>%
  dplyr::select(., -name)

# Put data and fitted parameters in same Df
countriesAll <- merge(countriesCovid, countriesFit, by="Id")

# calculate the fitted values for the countie dynamics
countriesAll <- countriesAll %>%
  dplyr::mutate(., t_bp = ceiling(10^time.threshold),
   cumI_fit = ifelse(t<10^time.threshold,
                     10^(intercept + log10(t)*slopeI),
                     10^(intercept - log10(t_bp)*slopeF +
                         log10(t_bp) * slopeI + log10(t) * slopeF)),
                Country = as.character(Country)) %>%
  as_tibble(.)


# select examples to show
increasingRatio <- dplyr::filter(countriesFit, slopeRatio > 1)
slopeOrder <- base::order(increasingRatio$slopeRatio)
nCountries <- length(slopeOrder)
#ind <- round(c(seq(1, nCountries, nCountries/(nExamples-1)), nCountries))
ind <- c(seq(1, nExamples/2, 1), seq(nCountries-nExamples/2+1, nCountries, 1))
examplesIndices <- slopeOrder[ind]
examplesId <- as.character(increasingRatio$Country[examplesIndices])

#plot
plotsCountries <- dplyr::filter(countriesAll, Country %in% examplesId) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#4020ab") +
  facet_wrap(~Country, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("time (days)") +
  #ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("Countries and regions") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")


plotAll <- ggarrange(plotsCounties, plotsCountries, nrow = 2)

saveRDS(list(countiesDyn = plotsCounties, countriesDyn = plotsCountries),
        "./generated_data/data_dynamics_plots.Rds")

ggsave(plotAll, file=fileAllDyns, width=5.5, height= 9, units = "in")
#ggsave(plot_condados_paises, file="fig_2_trayectorias.png",width=5.5,height=4)



#################
# Make plot S3 with the dynamics of all data
#################

countiesId <- as.character(countiesFit$Id)
nCounties <- length(countiesId)
plotCols <- 8
chunkRows <- 10
chunkPlots <- plotCols * chunkRows
nChunks <- ceiling(nCounties / chunkPlots)

for (nc in 1:nChunks) {
  indI <- (nc-1)*chunkPlots+1
  if (nc < nChunks) {
    indF <- nc*chunkPlots
  } else {
    indF <- length(countiesId)
    chunkRows <- ceiling((indF-indI+1)/plotCols)
  }
  chunkId <- countiesId[indI:indF]

  plotsCountiesSupp <- dplyr::filter(countiesAll, Id %in% chunkId) %>%
    ggplot(., aes(x = t, y = cases)) +
    geom_point(color = "#000080") +
    facet_wrap(~county.full, scales = "free", ncol = 8) +
    scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
    #scale_y_continuous(limits = c(10, NA), trans = "log10") +
    scale_y_continuous(trans = "log10") +
    xlab("time (days)") +
    #ylab("Cumulative infected") +
    ylab(element_blank()) +
    ggtitle("U.S. counties") +
    theme_classic() +
    theme(strip.background = element_rect(linetype = 0),
          #strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    geom_line(aes(x=t,y=cumI_fit), color="sky blue")

  fileName <- paste("./plots/S3_", as.character(nc), "all_counties_dyn.png", sep = "")
  ggsave(fileName, plotsCountiesSupp, width = 12,
         height = chunkRows*1.5, units = "in", limitsize = FALSE)
}

#plot
countriesId <- as.character(countriesFit$Country)
nCountries <- length(countriesId)
plotCols <- 8
chunkRows <- 10
chunkPlots <- plotCols * chunkRows
nChunks <- ceiling(nCountries / chunkPlots)

#fix names

for (nc in 1:nChunks) {
  indI <- (nc-1)*chunkPlots+1
  if (nc < nChunks) {
    indF <- nc*chunkPlots
  } else {
    indF <- length(countriesId)
    chunkRows <- ceiling((indF-indI+1)/plotCols)
  }
  chunkId <- countriesId[indI:indF]

  plotsCountriesSupp <- dplyr::filter(countriesAll, Country %in% chunkId) %>%
    dplyr::mutate(country = ifelse(Country_Region %in% "Korea, South","South Korea",
                                   Country_Region)) %>% 
    ggplot(., aes(x = t, y = cumI)) +
    geom_point(color = "#4020ab") +
    facet_wrap(~Country, scales = "free", ncol = plotCols) +
    scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
    #scale_y_continuous(limits = c(10, NA), trans = "log10") +
    scale_y_continuous(trans = "log10") +
    xlab("time (days)") +
    #ylab("Cumulative infected") +
    ylab(element_blank()) +
    ggtitle("Countries and regions") +
    theme_classic() +
    theme(strip.background = element_rect(linetype = 0),
          #strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size=10)) +
    geom_line(aes(x=t,y=cumI_fit),color="sky blue")

  fileName <- paste("./plots/S3_", as.character(nc), "all_countries_dyn.png", sep = "")
  ggsave(fileName, plotsCountriesSupp, width = 12,
         height = chunkRows*1.5, units = "in", limitsize = FALSE)
}

