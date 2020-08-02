#levanta datos para condados y países
require(stringi)
require(lubridate)
require(dplyr)
condPop = read.csv('counties.11.7.csv') %>% select(FIPS,Admin2,Province_State,Population)
colnames(condPop)[1:3]=c("fips","county","state")


casosCondados = read.csv('time_series_covid19_confirmed_US.csv') 
View(casosCondados)


us_counties = casosCondados %>% filter(! (Admin2 %in% "")) %>%  
  select(FIPS,Admin2,Province_State,starts_with("X")) %>%
  gather(key="date",value="cumI",4:174) %>%
  mutate(date=stri_replace(date,regex="^X",replacement = ""),date=mdy(date)) %>% 
  group_by(FIPS) %>% 
  mutate(neg=sum(diff(cumI)<0)) %>%
  filter(neg==0,cumI>0) 

colnames(us_counties)[1:3]=c("fips","county","state")



#casos.world.264 <- read.delim("casos.world.264.tab")
cases.world <- read.csv("world.data.tab.csv",sep = "\t")
head(cases.world)
#colnames(casos.world.264)
cases.world = cases.world %>% filter(!is.na(code3))
tail(cases.world$iso3)
tail(casos.world.264$iso3)
