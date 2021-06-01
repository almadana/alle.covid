source("1_Allee_phase_space.R")

I50=2*maxDetected
# with vaccination
# calculate Rt for given parameters
calculate_Rt_vacc <- function(maxDetected, maxLinks, Infected, pCall, NDailyCalls,
                              pInfection, Npop,pv){
  # Expected number of detected individuals
  pDetected <- maxDetected / (I50 + Infected)
  NDetected <- pDetected * Infected
  # Daily probability that a detected individual is contacted
  pContact <- 1 - (1-pCall)^(NDailyCalls/NDetected)
  # Number of expected social links for detected individuals
  detectedLinks <- expected_links(maxLinks, pContact)
  # Probability that a contact of an infected is susceptible
  pSusceptible<- (Npop - Infected - 1) / (Npop - 1)
  # Rt calculation
  Rt <- pv*pSusceptible * pInfection * (pDetected*detectedLinks + (1-pDetected)*maxLinks)
  return(Rt)
}

infected <- c(1:round(Npop*propInfectedMax))

# calculates Rt for a vector of values of Infected
calculate_Rt_range_vacc <- function(maxDetected, maxLinks, Infected, pCall,
                                    NDailyCalls, pInfection, Npop,vp){
  RtVec <- vector()
  for (i in c(1:length(Infected))) {
    RtVec[i] <- calculate_Rt_vacc(maxDetected, maxLinks, Infected[i], pCall,
                                  NDailyCalls, pInfection, Npop,vp)
  }
  return(RtVec)
}

vp=1
vp=c(0.3,0.6,1)
a=tibble(lapply(vp,function(x) calculate_Rt_range_vacc(maxDetected,maxLinks,infected,pCall,NDailyCalls,pInfection,Npop,x)))
a$vp = vp
a=a %>% unnest()
a$Npop = rep(1:Npop,3)
colnames(a)=c("Rt","fnv","Prop.I")

ggplot(a,aes(x=Prop.I,y=Rt,col=as.character(fnv)))+geom_point()+
  labs(x="Proportion infected",y="Rt",col="Proportion unvaccinated")