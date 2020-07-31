cuadrantes_long = cuadrantes %>% head(100) %>% 
  gather(key=t,value=cumI,starts_with("t.")) %>% 
  mutate(t=strotoi(stringr::str_replace(pattern="t.",replacement = "")))

plotsDynAllee <- dplyr::filter(cuadrantes, allee == TRUE & rep %in%sample4) %>%
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
