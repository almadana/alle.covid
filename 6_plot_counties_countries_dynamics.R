library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

fileAllDyns <- "./plots/data_dynamics.pdf"

fileSuppDyns <- "./plots/all_counties_dyn.pdf"

######################################
####### Plot counties dynamics #######
######################################

# load data
countiesFit <- readRDS("./generated_data/countiesFit.RDS") %>%
  dplyr::mutate(., slopeRatio = slopeF / slopeI) %>%
  as_tibble(.)

fittedCounties <- as.character(countiesFit$Id)

countiesCovid <- read.csv("./raw_data/us-counties.csv", colClasses = c("fips" = "character")) %>%
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
                         log10(t_bp-1) * slopeI + log10(t) * slopeF))) 

# get counties with the least and most pronounced breakpoints
nExamples <- 4 # cuantos de cada
increasingRatio <- dplyr::filter(countiesFit, slopeRatio > 1)
slopeOrder <- base::order(increasingRatio$slopeRatio)
nCounties <- length(slopeOrder)
examplesIndices <- slopeOrder[c(1:nExamples, (nCounties-nExamples+1):nCounties)]
examplesId <- as.character(countiesFit$Id[examplesIndices])

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
  dplyr::mutate(., slopeRatio = slopeF / slopeI) %>%
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
                         log10(t_bp-1) * slopeI + log10(t) * slopeF))) %>%
  as_tibble(.)


# select examples to show
nExamples <- 4 # cuantos de cada
increasingRatio <- dplyr::filter(countriesFit, slopeRatio > 1)
slopeOrder <- base::order(increasingRatio$slopeRatio)
nCountries <- length(slopeOrder)
examplesIndices <- slopeOrder[c(1:nExamples, (nCountries-nExamples+1):nCountries)]
examplesId <- as.character(countriesFit$Id[examplesIndices])

#plot
plotsCountries <- dplyr::filter(countriesAll, Id %in% examplesId) %>%
  dplyr::mutate(country = ifelse(Country_Region %in% "Korea, South","South Korea",
                                 Country_Region)) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#4020ab") +
  facet_wrap(~Country_Region, scales = "free", ncol = 4) +
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

  fileName <- paste("./plots/S3_", as.character(nc), "all_counties_dyn.pdf", sep = "")
  ggsave(fileName, plotsCountiesSupp, width = 22,
         height = chunkRows*2.5, units = "in", limitsize = FALSE)
}

#plot
countriesId <- as.character(countriesFit$Id)
nCountries <- length(countriesId)
plotCols <- 8
chunkRows <- 10
chunkPlots <- plotCols * chunkRows
nChunks <- ceiling(nCountries / chunkPlots)

for (nc in 1:nChunks) {
  indI <- (nc-1)*chunkPlots+1
  if (nc < nChunks) {
    indF <- nc*chunkPlots
  } else {
    indF <- length(countriesId)
    chunkRows <- ceiling((indF-indI+1)/plotCols)
  }
  chunkId <- countriesId[indI:indF]

  plotsCountriesSupp <- dplyr::filter(countriesAll, Id %in% chunkId) %>%
    dplyr::mutate(country = ifelse(Country_Region %in% "Korea, South","South Korea",
                                   Country_Region)) %>% 
    ggplot(., aes(x = t, y = cumI)) +
    geom_point(color = "#4020ab") +
    facet_wrap(~Country_Region, scales = "free", ncol = plotCols) +
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

  fileName <- paste("./plots/S3_", as.character(nc), "all_countries_dyn.pdf", sep = "")
  ggsave(fileName, plotsCountriesSupp, width = 20,
         height = chunkRows*2.2, units = "in", limitsize = FALSE)
}

