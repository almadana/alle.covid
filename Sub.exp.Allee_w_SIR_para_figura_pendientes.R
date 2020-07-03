###################################################################
##### Dynamic

dyn.sub.exp.Allee.SIR<-function(I0, bet=1.4, gam=5, p, t, Allee, migration, size, I50, K){
  out<-matrix(0, nrow = t, ncol=5)
  colnames(out)<-c("S","I","R", "New.I", "Imported")
  out[1,]<-c(K,I0,0,0,0)
  N<-K+I0
  for(i in 2:t){
    if(out[i-1,2]<=0)out[i-1,2]<-out[i,2]<-0
    if(out[i-1,1]<=0)out[i-1,1]<-out[i,1]<-0
        if(Allee==T){
          new.I<-(bet*(out[i-1,1]*(out[i-1,2])^p)/N)*(out[i-1,2]/(out[i-1,2]+I50))
          if(out[i-1,2]>0){
              rec.I<-(1/(gam*(out[i-1,2]/(out[i-1,2]+I50))))*out[i-1,2]
              } else {rec.I<-0}
          }
          if(Allee==F){
            new.I<-(bet*(out[i-1,1]*(out[i-1,2])^p)/N)
            rec.I<-(1/gam)*out[i-1,2]
        }
        rec.I<-rnbinom(n = 1, size = size, mu = rec.I)
        rec.I<-min(out[i-1,2],rec.I, na.rm=T)
        imported<-rnbinom(n = 1, size = size, mu = migration)
#        new.I<-rnbinom(n = 1, size = size, mu = new.I)
        new.I<-rpois(n = 1, lambda  = new.I)
        new.I<-min(new.I,out[i-1,2])        # ensure no more new infections than susceptible
        out[i,]<-c(out[i-1,1]-new.I, out[i-1,2]+new.I-rec.I+imported,out[i-1,3]+rec.I, new.I, imported)
    }
  out
}


##############################################
#############################################

p_se=1 # sin subexponencial
#p_se=0.8 # subexponencial
allee=T

par(mfrow=c(10,10), mar=c(1.5,1.5,0.5,0.5 ))
#for(i in 1:200){
generador <- function(p_se,allee) {
  dyn.sub.exp.Allee.SIR(I0=10, bet=0.8, gam=3 , p=p_se, t=100,  migration=1, size=0.2, I50=10, K=100000, Allee = allee )->al
  imp<-cumsum(al[,5])  # cumulative number of imported cases
  al<-cumsum(al[,4])   # cumulative number of cases
  I.min<-which(al>9)
  al<-al[I.min]
  I.max<-which(al<0.8*1000)
  al<-al[I.max]
  print(length(al))
  al<-al[which(al>10 & al<2500)]
  imp<-imp[which(al>10 & al<2500)]
  #al<-al[which(al>10 )]
  #imp<-imp[which(al>10 )]
  if(length(al)>10){
    plot(log10(al)~log10(I(1:length(al))), bty="l", pch=19, col="darkgreen", main="")
  y<-log10(al)
  z<-log10(imp)
  x<-log10(1:length(y))
#  points(z~log10(I(1:length(al))), bty="l", pch=19, col="gold", main="", cex=.8)
  #ff0<-lm(y~x)
  ff1<-segmented(lm(y~x), seg.Z = ~x, npsi = 1)
  #aic<-AIC(ff0,ffp, ff1)[,2]
  #wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
  cfs = summary(ff1)
  
  return(c(cfs$coefficients[,1],ifelse(is.null(cfs$psi),NA,cfs$psi[2])))
#  plot(ff1, add=T, col="red")
  }
}
# 
# par(mfrow=c(1,1), mar=c(1.5,1.5,1.5,1.5 ))
# 
# d=data.frame(n=1:5000)
# fits=sapply(d$n,FUN = generador)
# fnull=sapply(fits,FUN = is.null)
# sum(fnull)
# fits=fits[!fnull]
# 
# 
# cociente=sapply(fits, FUN = function(x) {
#   cfs = coef(x)
#   return(cfs[3] / cfs[2])
# })
# hist(log(cociente),20)
# 
# resta=sapply(fits, FUN = function(x) {
#   cfs = coef(x)
#   return(cfs[3] - cfs[2])
# })
# hist(resta,20)
# 

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(segmented)
nRepeats = 2000
alleeCoefs = c(F,T)
#subExpCoefs = c(0.8,1) #con y sin subexp
subExpCoefs = 0.8
nRows=length(alleeCoefs)*length(subExpCoefs)

cuadrantes = tibble(n=seq(1,nRows*nRepeats),
                  allee=rep(alleeCoefs,nRows/length(alleeCoefs)*nRepeats),
                    allee.f = ifelse(allee,"with Allee effect","without Allee effect"),
                    subexp = rep(subExpCoefs,each=nRows/length(subExpCoefs)*nRepeats),subexp.f = ifelse(subexp<1,"with subexp. growth","without subexp. growth"))

cuadrantes = cuadrantes %>% rowwise() %>% 
  # mutate(modelo=list(generador(subexp,allee)),fit=list(coef(modelo,include.psi=T)))
  mutate(fit=list(generador(subexp,allee)))

cuadrantes = 
cuadrantes %>% 
  unnest(fit)  %>% 
  group_by(n) %>% 
  mutate(coef=paste0("coef.",1:n())) %>% 
  spread(key=coef,value  = fit) %>% 
  rowwise() %>% 
  mutate(cociente = coef.3/coef.2,resta=coef.3 - coef.2,
         angle = atan(abs((coef.3 - coef.2)/(1-coef.3*coef.2))),
         infec.Corte=ifelse(is.null(coef.5),NA,coef.1 + coef.2*coef.5)) 


#cuadrantes = cuadrantes %>% rowwise() %>% mutate(angle = atan(abs((coef.3 - coef.2)/(1-coef.3*coef.2)))  )
#cuadrantes = cuadrantes %>%ungroup() %>%  rowwise() %>% mutate(infec.Corte = coef.1 + coef.2*coef.5  )

#cuadrantes$coef.1 + cuadrantes$coef.2*cuadrantes$coef.5

cuadrantes %>% group_by(allee.f,subexp.f) %>% 
  summarize_at(c("cociente","resta","angle"),.funs = list("median"=median,"iqr"=IQR),na.rm=T)

#filtrar casos con muy poco slope inicial y mucho final - cociente>5000
cuadrantes=cuadrantes %>%  filter(cociente<5000,cociente>0)

cuadrantes %>% ggplot(aes(x=log(cociente),fill=interaction(allee.f,subexp.f))) + geom_histogram(alpha=.3) + facet_wrap(~subexp.f*allee.f)
cuadrantes %>% ggplot(aes(x=log(cociente),fill=allee)) + geom_histogram() 
cuadrantes %>% ggplot(aes(x=resta,fill=allee.f)) + geom_histogram() + facet_grid(rows = "allee.f~subexp.f")

cuadrantes %>% ggplot(aes(x=log(cociente),fill=allee.f)) + geom_histogram() + facet_grid(rows = "allee.f~subexp.f",)

cuadrantes %>% ggplot(aes(y=log(cociente),x=log(1+coef.2),col=allee.f)) + geom_point(alpha=.2)+
  facet_wrap(~subexp.f)+
  labs(x="log(1+initial slope)",y="log(second slope/first slope)")

cuadrantes %>% ggplot(aes(y=log(1+coef.3),x=log(1+coef.2),col=allee.f)) + geom_point(alpha=.2)+
  facet_wrap(~subexp.f)+
  labs(x="log(first slope)",y="log(second slope)")

cuadrantes$allee.f = factor(cuadrantes$allee.f)

# cuadrantes %>% mutate_at(c("coef.2","coef.3","cociente"),.funs = list("log"=function(x) log(1+x))) %>% 
# ggscatterhist(x="coef.2_log",y="cociente_log",
#   color="allee.f",alpha=0.6,size=3,
#   margin.plot="boxplot",
#   ggtheme=theme_bw(),
#   xlab="log(1+initial slope)",
#   ylab="log(1+final slope / initial slope)"
# ) 

cuadrantes %>% mutate_at(c("coef.2","coef.3","cociente"),.funs = list("log10"=function(x) log10(1+x))) %>% 
ggscatterhist(x="coef.2_log10",y="coef.3_log10",
              color="allee.f",alpha=0.6,size=3,
              margin.plot="boxplot",
              ggtheme=theme_bw(),
              xlab="log10(1+initial slope)",
              ylab="log10(1+final slope)"
)

cuadrantes %>% mutate_at(c("coef.2","coef.3","cociente"),.funs = list("log"=function(x) log(1+x))) %>% 
  ggscatterhist(x="coef.2",y="angle",
                color="allee.f",alpha=0.6,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="initial slope",
                ylab="angle between slopes"
  )

cuadrantes %>% 
  ggscatterhist(x="infec.Corte",y="angle",
                color="allee.f",alpha=0.6,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(infected at threshold)",
                ylab="angle between slopes"
  )

cuadrantes %>% 
  ggscatterhist(x="coef.5",y="angle",
                color="allee.f",alpha=0.6,size=3,
                margin.plot="boxplot",
                ggtheme=theme_bw(),
                xlab="log10(time to threshold)",
                ylab="angle between slopes"
  )





cuadrantes %>% ggplot(aes(y=coef.3,x=coef.2,col=allee.f)) + geom_point(alpha=.2)+
  facet_wrap(~subexp.f)+
  labs(x="initial slope)",y="final slope")


cuadrantes %>% ggplot(aes(y=log(cociente),col=interaction(allee,subexp.f),x=resta)) + geom_point() +
  facet_wrap(~subexp.f)


f###################################################################
##### plot with/wihout Allee and subexponential

plot.alle.migra.su.exp<-function(I0, bet, gam , p, t,  migration, size, I50, K, it){
  require(segmented)
  par(mfrow=c(2,2), mar=c(4,4,4,4))
  dyn.median<-list()
  for(pp in c(1,p)){
    #    for(m in c(0,migration)){
    for(allee in c(F,T)){
      out.t<-NULL
      for(i in 1:it){
        dyn.sub.exp.Allee.SIR(I0 = I0, bet = bet, gam = gam, p = pp, t = t, Allee = allee, migration = migration, size = size, I50 = I50, K = K)->al
        #return(al)       
        al<-cumsum(al[,4])
        out.t<-rbind(out.t, al)
      }
      out<- apply(out.t, 2, quantile, c(0.025,0.5,0.975))
      #    ii<-which(diff(al)==max(diff(al)))
      Total.Infected<-out[2,]
      
      if(sum(Total.Infected)==0)return("No epidemic")      
      if(allee ==T & pp<1 )          title<-"With Allee & With Sub.exp"
      if(allee ==T & pp==1 )         title<-"With Allee & Without Sub.exp "
      if(allee ==F & pp==1 )         title<-"Without Allee & Without Sub.exp "
      if(allee ==F & pp<1 )          title<-"Without Allee & With Sub.exp"
      if(allee==F){  
        Total.Infected<-Total.Infected[which(Total.Infected>1 & Total.Infected<K*0.25)]
        #         dyn.median[length(dyn.median)+1]<-c(allee, pp,Total.Infected)
        plot(log10(Total.Infected)~log10(I(1:length(Total.Infected))), bty="l", pch=19, 
             col="navy", main=title, ylim=c(0.95,max(log10(Total.Infected))),
             xlab="Log10(day)", ylab="Log10(total Infected)", xlim=c(0,log10(t)))
        #      points(log10(out[1,which(Total.Infected>9)])~log10(I(1:length(Total.Infected))), col="navy", lty=2, type="l")
        #      points(log10(out[3,which(Total.Infected>9)])~log10(I(1:length(Total.Infected))), col="navy", lty=2, type="l")
      }
      if(allee==T){
        Total.Infected<-Total.Infected[which(Total.Infected>1 & Total.Infected<K*0.25)]
        #        dyn.median<-rbind(dyn.median, c(allee, pp,Total.Infected))
        plot(log10(Total.Infected)~log10(I(1:length(Total.Infected))), bty="l", pch=19, col="darkred", main=title)
        # points(log10(out[1,I.min])~log10(I(1:length(Total.Infected))), col="darkred", lty=2, type="l")
        #  points(log10(out[3,I.min])~log10(I(1:length(Total.Infected))), col="darkred", lty=2, type="l")
        y<-log10(Total.Infected); x<-log10(1:length(y))
        DI<-log10(diff(Total.Infected))
        ff0<-lm(y~x)
        ff1<-segmented(lm(y~x), seg.Z = ~x, npsi = 1)
        ff3<-lm(y~x+I(x^2)+I(x^3))
        #        #Vivouda et al. 2016 Epidemics 15: 27–37 (DI=r*Ct^p)
        aic<-AIC(ff0,ff1, ff3)[,2]
        wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
        print(wi)
        fit<-summary(ff1)
        plot(ff1, add=T, col="red")
      }
    }
  }
  #colnames(dyn.median)<-c("Allee", "sub.exp",1:t)  
  dyn.median
}




#I0 = 0, bet = 1.2, gam = 5.5, p = 0.8, t = 100, Allee = T, migration = 3, size = 0.75, I50 = 20, K = 10000

#dd<-plot.alle.migra.su.exp(I0=0, bet=1.2, gam=5.5 , p=.8, t=200,  migration=1, size=0.5, I50=20, K=10000, it=200)
#dd<-plot.alle.migra.su.exp(I0=0, bet=0.8, gam=3.5 , p=.8, t=100,  migration=1, size=0.25, I50=20, K=10000, it=500)
dd<-plot.alle.migra.su.exp(I0=0, bet=0.8, gam=3 , p=.8, t=200,  migration=1, size=0.2, I50=20, K=100000, it=500)


