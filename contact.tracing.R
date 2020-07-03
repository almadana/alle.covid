# L.max, is the maximum number of effictive contacts that an infected can achieve (e.g. having enough time to effectively transmit the virus)
# Infected, is the number of infected people in the population
# B.link, is the proibablity of transmision through an effective link. It is determined by behaviour as mask use, persons distances, etc
# N.tree.day, is the number of contacts tree that can be majke in a day
# n.by.tree, is the average number of individuals in a contact tree
# day.to.report, average number of days from the INFECTION to the report (detection by syntoms, contact tracing, or test)

# P.find, is the probability of find the exposed person of the tree (in a phone call or app for cell tracking)
# N.population is the total population size
# alpha is the fraction of asyntomatics

##LINEAS PARA CORRER: a<-space.Rt.inf(L.max = 12, P.find = 0.2, n.by.tree = 10, B.link = 0.3, day.to.report = 7, N.population = 2000 )
#b<-ifelse(a>1,1,0)
#image(t(b[1:1600,1:1000]), col=c("royalblue2", "firebrick3")) 


#p, subexponential parameter

R.time<-function(L.max,Infected, P.find, N.call.day, N.tree.day, n.by.tree, B.link, day.to.report, N.population){

  Max.N<-N.tree.day*n.by.tree                                 # Maximum number of individulas that can be identfied in the conact tree
  detected<-Infected*Max.N/(Max.N+Infected)                   # detected cases is the expected Binomial with N: number of infected, p: 
                                                              # Probability of detection (is a decreasing function of infected cases (Hill function with exponent 1)) 
  
  P.effective.call<-1-(1-P.find)^(N.call.day/detected)        # is the probability of contact and alert the exposed. It involve the probability of fin th person
                                                              # in the contact (e.g. phone or app)

  Exp.day<- 1/P.effective.call                               # Expected day from a geometric distribution (day at which stop contacting)
  Links<-(L.max*(Exp.day^4)/((Exp.day^4)+7^4))               # Number of effective contacts. At day 7 half of the maximum effective contacts (links) take place. Until day 4 no effective links
  f.sus<-(N.population-Infected)/N.population        # OJO CAMBIR ENmodelo dinamico POR (S-I-R)/S =f.sus probability that he linked indiviual is susceptible
  Rt<-((f.sus)*((detected/Infected)*Links*B.link +((Infected-detected)/Infected)*L.max*B.link))  #cambiar por estimacion en base a dias de infectado....  en modelo dinamico Fsus
  Rt
}

R.time.2.0<-function(L.max,Infected, P.find, N.call.day, N.tree.day, n.by.tree, B.link, day.to.report, N.population){
  detected<-Infected*Max.N/(Max.N+Infected)                   # detected cases is the expected Binomial with N: number of infected, p: 
  # Probability of detection (is a decreasing function of infected cases (Hill function with exponent 1)) 
  
  P.detected<-P.find*min(1.)*(1-(1-P.find)^(P.find*Infected/detected))        # is the probability of contact and alert the exposed. It involve the probability of fin th person
  # in the contact (e.g. phone or app)
  
  Exp.day<- 1/P.detected                               # Expected day from a geometric distribution (day at which stop contacting)
  Links<-(L.max*(Exp.day^4)/((Exp.day^4)+7^4))               # Number of effective contacts. At day 7 half of the maximum effective contacts (links) take place. Until day 4 no effective links
  f.sus<-(N.population-Infected)/N.population        # OJO CAMBIR ENmodelo dinamico POR (S-I-R)/S =f.sus probability that he linked indiviual is susceptible
  Rt<-((f.sus)*((detected/Infected)*Links*B.link +((Infected-detected)/Infected)*L.max*B.link))  #cambiar por estimacion en base a dias de infectado....  en modelo dinamico Fsus
  Rt
  
  
}


##############################
# Space Rt~ Infected + strategy

### Number of calls
space.Rt.inf<-function(L.max, P.find,  n.by.tree, B.link, day.to.report, N.population){
  out<-matrix(NA,ncol=800, nrow=1500)
  for (i in 1:1500){
print(i)    
    for(j in 1:800){

    rr<-R.time(L.max=L.max,Infected=i, P.find=P.find, N.call.day=j, N.tree.day=max(round(j*0.2),1), n.by.tree=n.by.tree, B.link=B.link, day.to.report=day.to.report, N.population=N.population)
    out[i,j]<-rr
    }
  }
out  
}

################

space.Rt.inf.L.max<-function(P.find,  n.by.tree, B.link, day.to.report, N.population, N.call.day, N.tree.day){
  out<-matrix(NA,ncol=100, nrow=2000)
  for (i in 1:2000){
    print(i)    
    for(j in 1:100){
      
      rr<-R.time(L.max=j,Infected=i, P.find=P.find, N.call.day=N.call.day, N.tree.day=N.tree.day, n.by.tree=n.by.tree, B.link=B.link, day.to.report=day.to.report, N.population=N.population)
      out[i,j]<-rr
    }
  }
  out  
}

#############
space.Rt.inf.Bt<-function(L.max,P.find,  n.by.tree, day.to.report, N.population, N.call.day, N.tree.day){
  out<-matrix(NA,ncol=300, nrow=2000)
  Bl<-seq(0,1,,300)
  for (i in 1:2000){
    print(i)    
    for(j in 1:300){
      rr<-R.time(L.max=L.max,Infected=i, P.find=P.find,
                 N.call.day=N.call.day, N.tree.day=N.tree.day, 
                 n.by.tree=n.by.tree, B.link=Bl[j], day.to.report=day.to.report,
                 N.population=N.population)
      out[i,j]<-rr
    }
  }
  out  
}






