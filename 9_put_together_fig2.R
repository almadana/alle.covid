library(ggplot2)
library(ggpubr)
source("./7_scatter_counties.R")
#source("./8_scatter_countries.R")
source("./3_plot_Allee_dyn.R")
# load and organize the individual plots into figure 2
#dynamicsPlots <- readRDS("./generated_data/data_dynamics_plots.Rds")

fig2 <- ggarrange(dynsPlots,
                  scatterCountiesPlot,
                  # dynamicsPlots$countiesDyn, 
                  # scatterCountriesPlot,
                  # dynamicsPlots$countriesDyn,
                  widths = c(6,4),
                  labels = "AUTO")

ggsave(fig2, file="./plots/fig2.png", width=14, height=6, units = "in")
ggsave(fig2, file="./plots/fig2.pdf", width=14, height=6, units = "in")

