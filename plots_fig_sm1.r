# FIGURA CONDADOS Y PAÍSES SUPP MAT

# ----- FIGURE S1 - COUNTIES -------

# cuando quise poner una función de ploteo no me anduvo :P


countyList = county_data_fit %>% arrange(county.full) %>% group_by(county.full) %>% slice(1) %>%  pull(Id)



#n_per_fig = round(length(countyList)/3)
#al final 80 80 y 76 anda bien


#plotsCondados_SM <- county_data_fit %>% filter(Id %in% countyList[1:n_per_fig]) %>% 
plotsCondados_SM <- county_data_fit %>% filter(Id %in% countyList[1:80]) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#000080",size=.6) +
  facet_wrap(~county.full, scales = "free", ncol = 8) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  ggtitle("U.S. counties") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size=8),
        text = element_text(size=5)) +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

ggsave(plotsCondados_SM,file="sm1_1.pdf",width=6,height = 8)
ggsave(plotsCondados_SM,file="sm1_1.png",width=6,height = 8)


plotsCondados_SM <- county_data_fit %>% filter(Id %in% countyList[81:160]) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#000080",size=.6) +
  facet_wrap(~county.full, scales = "free", ncol = 8) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  #xlab("time (days)") +
  #ylab("Cumulative infected") +
  #ylab(element_blank()) +
  ggtitle("U.S. counties") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=5)) +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")


ggsave(plotsCondados_SM,file="sm1_2.pdf",width=6,height = 8)
ggsave(plotsCondados_SM,file="sm1_2.png",width=6,height = 8)


plotsCondados_SM <- county_data_fit %>% filter(Id %in% countyList[161:228]) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#000080",size=.6) +
  facet_wrap(~county.full, scales = "free", ncol = 8) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  #xlab("time (days)") +
  #ylab("Cumulative infected") +
  #ylab(element_blank()) +
  ggtitle("U.S. counties") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=5)) +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")


ggsave(plotsCondados_SM,file="sm1_3.pdf",width=6,height = 6.7)
ggsave(plotsCondados_SM,file="sm1_3.png",width=6,height = 8)



## --  FIGURE S2 - COUNTRIES ---------



countryList = country_data_fit %>% arrange(country) %>% group_by(country) %>% slice(1) %>%  pull(Id)



#n_per_fig = round(length(countryList)/3)
#al final 80 80 y 76 anda bien


plotsPaises_SM <- country_data_fit %>% filter(Id %in% countryList[1:80]) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#4020ab",size=.6) +
  facet_wrap(~country, scales = "free", ncol = 8) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  ggtitle("Countries & regions") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size=8),
        text = element_text(size=5)) +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

ggsave(plotsPaises_SM,file="sm2_1.pdf",width=6,height = 8)
ggsave(plotsPaises_SM,file="sm2_1.png",width=6,height = 8)


plotsPaises_SM <- country_data_fit %>% filter(Id %in% countryList[81:113]) %>% 
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#4020ab",size=.6) +
  facet_wrap(~country, scales = "free", ncol = 8) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  #xlab("time (days)") +
  #ylab("Cumulative infected") +
  #ylab(element_blank()) +
  ggtitle("Countries & regions") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=5)) +
  xlab(element_text("time (days)",size=10)) +
  ylab(element_text("Cumulative infected",size=10)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")


ggsave(plotsPaises_SM,file="sm2_2.pdf",width=6,height = 4)
ggsave(plotsPaises_SM,file="sm2_2.png",width=6,height = 4)



