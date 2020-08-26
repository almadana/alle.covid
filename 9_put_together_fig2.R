library(ggplot2)
library(ggpubr)
source("./7_scatter_counties.R")
source("./8_scatter_countries.R")

# load and organize the individual plots into figure 2
dynamicsPlots <- readRDS("./generated_data/data_dynamics_plots.Rds")
scatterCounties <- readRDS("./generated_data/counties_scatter_plot.RDS")
scatterCountries <- readRDS("./generated_data/countries_scatter_plot.RDS")

fig2 <- ggarrange(dynamicsPlots$countiesDyn,
          scatterCountiesPlot, 
          dynamicsPlots$countriesDyn,
          scatterCountriesPlot,
          widths = c(4,3),
          labels = "auto")

ggsave(fig2, file="./plots/fig2.pdf", width=14, height=9, units = "in")
#ggsave(fig2, file="fig2.png",width=14,height=9)

