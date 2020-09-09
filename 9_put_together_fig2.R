library(ggplot2)
library(ggpubr)
source("./7_scatter_counties.R")
source("./8_scatter_countries.R")

# load and organize the individual plots into figure 2
dynamicsPlots <- readRDS("./generated_data/data_dynamics_plots.Rds")

fig2 <- ggarrange(scatterCountiesPlot,
                  dynamicsPlots$countiesDyn, 
                  scatterCountriesPlot,
                  dynamicsPlots$countriesDyn,
                  widths = c(4,4),
                  labels = "auto")

ggsave(fig2, file="./plots/fig2.png", width=10, height=6.5, units = "in")
ggsave(fig2, file="./plots/fig2.pdf", width=10, height=6.5, units = "in")

