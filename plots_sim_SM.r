nPlotDynSM=40
sampleRepetitionsSM <- sample(unique(simulations$rep), nPlotDynSM)


plotsDynAlleeSM <- dplyr::filter(simulations_fit, allee == TRUE & rep %in% sampleRepetitionsSM) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#b33018",size=.6) +
  #  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  #ylab("Cumulative infected") +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  ggtitle("with Allee effect") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
  plot.title = element_text(hjust = 0.5,size=8),
  text = element_text(size=5)) +
  geom_line(aes(x=t,y=cumI_fit),color="pink")

plotsDynAlleeSM

#plotsDynNonAllee <- dplyr::filter(simulations_fit, allee == FALSE & rep %in% sampleRepetitions) %>%
plotsDynNonAlleeSM <- dplyr::filter(simulations_fit, allee == FALSE & rep %in% sampleRepetitionsSM) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#14b74b",size=.6) +
  #  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  #xlab("time (days)",size=10) +
  #ylab("Cumulative infected") +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  ggtitle("without Allee effect") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size=8),
        text = element_text(size=5)) +
  geom_line(aes(x=t,y=cumI_fit),color="green")

fig_sim_SM = ggarrange(plotsDynAlleeSM,plotsDynNonAlleeSM,ncol = 2,nrow=1)

ggsave("simDynamics_SM.png", fig_sim_SM, width = 6, height = 8)
ggsave("simDynamics_SM.pdf", fig_sim_SM, width = 6, height = 8)
  