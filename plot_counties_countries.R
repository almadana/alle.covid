#figura para plotear condaos y países

max(bb$Id)

head(us.counties)
head(nn)

bb$Id

#abreviaciones de estados
state_abbrev = read.csv('county_abbreviations.csv')

county_data =   
  us_counties %>% group_by(fips) %>% 
  mutate(t=1:n(),peak = which.max(diff(cumI))[1]) %>% 
  filter(t>13,t<=peak) %>% 
  mutate(t=t-13)
colnames(county_data)[1] ="Id"

#agregar nombre extendido de condado
county_data = merge(county_data,state_abbrev,by="state",all.x = T)
county_data$county.full = paste(county_data$county,county_data$code,sep=", ")

#bb sale de plotCondados, que a su vez sale de segmented_condados (convierte el ajuste b en un tible bb, corrije segunda pendiente y agrega cociente para plotear)
county_data_fit = merge(county_data,bb,by="Id")


colnames(county_data_fit)[17]="breakPoint"

county_data_fit = 
  county_data_fit %>% mutate(t_bp = ceiling(10^breakPoint),cumI_fit = ifelse(t<10^breakPoint, 10^(intercept + log10(t)*initial.slope),10^(intercept - log10(t_bp)*slope.after.thresh +log10(t_bp-1)*initial.slope +log10(t)*slope.after.thresh ))) 

#los condados con más y con menos quiebre:
nExamples=5 # cuantos de cada
county_examples = 
  c(
    county_data_fit %>% filter(cociente>1)%>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_max(cociente,n=nExamples) %>% pull(Id)
  ,
  county_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
    slice_min(cociente,n=nExamples) %>% pull(Id)
  )

county4 = county_examples[order(county_examples)]
county4 = county4[-c(4,9)]


plotsCondados <- dplyr::filter(county_data_fit, Id %in% county4) %>%
  ggplot(., aes(x = t, y = cumI)) +
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
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

plotsCondados





country_data  = 
cases.world %>% mutate(Id=1:n())  %>% group_by(UID) %>%
  select(-X) %>% 
  select(Id,Province_State,Country_Region,starts_with("X")) %>% 
  pivot_longer(starts_with("X")) %>%
  filter(value>0) %>% 
  mutate(t=1:n(),peak = which.max(diff(value))[1]) %>% 
  filter(t>14,t<=peak) %>%
  mutate(t=t-14) %>%
  ungroup() %>% rowwise() %>% 
  mutate(Province_State = ifelse(Province_State %in% "",levels(Country_Region)[Country_Region],levels(Province_State)[Province_State])) 
colnames(country_data)[c(3,6)] = c("country","cumI")


#bb_p sale de plotPaises, que a su vez sale de segmented_paises (convierte el ajuste b en un tible bb_p, corrije segunda pendiente y agrega cociente para plotear)
country_data_fit = merge(country_data,bb_p,by="Id")
#View(country_data_fit)

colnames(country_data_fit)[14]="breakPoint"

country_data_fit = 
  country_data_fit %>% mutate(t_bp = floor(10^breakPoint),cumI_fit = ifelse(t<10^breakPoint, 10^(intercept + log10(t)*initial.slope),10^(intercept - log10(t_bp)*slope.after.thresh +log10(t_bp-1)*initial.slope +log10(t)*slope.after.thresh ))) 




#los condados con más y con menos quiebre:
nExamples=5 # cuantos de cada
country_examples = 
  c(
    country_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_max(cociente,n=nExamples) %>% pull(Id)
    ,
    country_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_min(cociente,n=nExamples) %>% pull(Id)
  )


country_examples
unique(country_data_fit$country[ country_data_fit$Id %in% country_examples])

country4 = country_examples[order(country_examples)]
country4 = country4[-c(4,8)]


#
plotsPaises <- dplyr::filter(country_data_fit, Id %in% country4) %>%
#plotsPaises <- dplyr::filter(country_data_fit, Id %in% 108) %>%
  mutate(country=ifelse(country %in% "Korea, South","South Korea",country)) %>% 
#plotsPaises <- dplyr::filter(country_data_fit, Id %in% c(166,39)) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#008080") +
  facet_wrap(~country, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("time (days)") +
  #ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("Countries and regions") +
  theme_classic() +
    theme(strip.background = element_rect(linetype = 0),
#    strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=t,y=cumI_fit),color="blue")

plotsPaises

plot_condados_paises = ggarrange(plotsCondados,plotsPaises,nrow = 2)
ggsave(plot_condados_paises,file="fig_2_trayectorias.pdf",width=5,height=4)
ggsave(plot_condados_paises,file="fig_2_trayectorias.png",width=5,height=4)
