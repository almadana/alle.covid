###################################################################
##### Dynamic

dyn.sub.exp.Allee.SIR<-function(I0, betamax=1.3, gammamax=5, p, t, Allee, migration, size, I50, K){
  out<-matrix(0, nrow = t, ncol=5)
  colnames(out)<-c("S","I","R", "New.I", "Imported")
  out[1,]<-c(K,I0,0,0,0)
  N<-K+I0
  for(i in 2:t){
    if(out[i-1,2]<=0)out[i-1,2]<-out[i,2]<-0
    if(out[i-1,1]<=0)out[i-1,1]<-out[i,1]<-0
        if(Allee==T){
          new.I<-(betamax*(out[i-1,1]*(out[i-1,2])^p)/N)*(out[i-1,2]/(out[i-1,2]+I50))
          if(out[i-1,2]>0){
              rec.I<-(1/(gammamax*(out[i-1,2]/(out[i-1,2]+I50))))*out[i-1,2]
              } else {rec.I<-0}
          }
          if(Allee==F){
            new.I<-(betamax*(out[i-1,1]*(out[i-1,2])^p)/N)
            rec.I<-(1/gammamax)*out[i-1,2]
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

par(mfrow=c(2,2), mar=c(1.5,1.5,0.5,0.5 ))
#for(i in 1:200){
generador <- function(p_se,allee,N,migration=1,s=0.2,betamax=1.3,gammamax=5) {
  rm(.Random.seed, envir=globalenv())
  dyn.sub.exp.Allee.SIR(I0=10, betamax=betamax, gammamax=gammamax , p=p_se, t=100,  migration, s, I50=10, K=N, Allee = allee )->al
  imp<-cumsum(al[,5])  # cumulative number of imported cases
  al<-cumsum(al[,4])   # cumulative number of cases
  I.min<-which(al>9)
  al<-al[I.min]
  I.max<-which(al<0.8*1000)
  al<-al[I.max]
  #print(length(al))
  al<-al[which(al>10 & al<2500)]
  imp<-imp[which(al>10 & al<2500)]
  #al<-al[which(al>10 )]
  #imp<-imp[which(al>10 )]
  if(length(al)>10){
#    plot(log10(al)~log10(I(1:length(al))), bty="l", pch=19, col="darkgreen", main="")
    #plot(log10(al)~(I(1:length(al))), bty="l", pch=19, col="darkgreen", main="")
    
  y<-log10(al)
  z<-log10(imp)
  x<-log10(1:length(y))
#  x<-(1:length(y))
#  points(z~log10(I(1:length(al))), bty="l", pch=19, col="gold", main="", cex=.8)
  ff0<-lm(y~x)
  ff1<-segmented(lm(y~x), seg.Z = ~x, npsi = 1)
  aic<-AIC(ff0, ff1)[,2]
  wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
  cfs = summary(ff1)
  coeffs = cfs$coefficients
  if (length(coeffs)<3) {
    coeffs = c(coeffs,NA,NA)
    wi=c(1,0)
  }
  #plot(ff1, add=T, col="red")
  id.psi<-which(x<cfs$psi[2])                                              # The following lines estimate the number of active cases at the break point
  if(length(id.psi)>9) I.active<-al[max(id.psi)]-al[max(id.psi)-9]  # If the time series before the break is larger than 10 days retain the total infected 10 days before...
  if(length(id.psi)<=9)I.active<-al[max(id.psi)]                        # If the time series is shoreter take the total number of infected before the break
  
  return(c(coeffs,ifelse(is.null(cfs$psi),NA,cfs$psi[2]),al,I.active,wi[2]))
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
library(ggforce)
nRepeats = 2000
alleeCoefs = c(F,T)
pob_condados = c(4.49,0.68)
pob_paises = c(6.66,1.07)


#N=round(10^rnorm(nRepeats,pob_paises[1],pob_paises[2]))  #valores para paises
N=round(10^rnorm(nRepeats,pob_condados[1],pob_condados[2]))  #valores para condados
#subExpCoefs = c(0.8,1) #con y sin subexp
#para ajustar a países y condados
#subExpCoefs = 0.7
nRows=length(alleeCoefs)
#para explorar
subExpCoefs = runif(nRepeats)*0.2 + 0.7 #uniforme entre 0.7 y 0.9
betamax = runif(nRepeats)*0.11 + 1.07
gammamax = runif(nRepeats)*3 + 4
cuadrantes = tibble(n=seq(1,nRows*nRepeats),
                  allee=rep(alleeCoefs,nRows/length(alleeCoefs)*nRepeats),
                    allee.f = ifelse(allee,"with Allee effect","without Allee effect"),
                    subexp = rep(subExpCoefs,each=nRows/length(subExpCoefs)*nRepeats),subexp.f = ifelse(subexp<1,"with subexp. growth","without subexp. growth"),
                    N=rep(N,each=nRows),
                  betamax=rep(betamax,each=nRows/length(subExpCoefs)*nRepeats),
                  gammamax=rep(gammamax,each=nRows/length(subExpCoefs)*nRepeats))

cuadrantes = cuadrantes %>% rowwise() %>% 
  # mutate(modelo=list(generador(subexp,allee)),fit=list(coef(modelo,include.psi=T)))
  mutate(fit=list(generador(subexp,allee,N,betamax,gammamax)))

cuadrantes = 
cuadrantes %>% 
  unnest(fit)  %>% 
  group_by(n) %>% 
  mutate(coef=c(paste0("coef.",1:5),paste0("t.",stringr::str_pad(1:(n()-7),3,pad="0")),"I.active","Weighted.evidence")) %>% #View()
  spread(key=coef,value  = fit) %>% 
  rowwise() %>%
  mutate(coef.4 = coef.2+coef.3,
         cociente = coef.4/coef.2,
         resta=coef.4 - coef.2,
         angle = atan(abs((coef.4 - coef.2)/(1-coef.4*coef.2))),
         infec.Corte=ifelse(is.null(coef.5),NA,coef.1 + coef.4*coef.5))


#cuadrantes = cuadrantes %>% rowwise() %>% mutate(angle = atan(abs((coef.3 - coef.2)/(1-coef.3*coef.2)))  )
#cuadrantes = cuadrantes %>%ungroup() %>%  rowwise() %>% mutate(infec.Corte = coef.1 + coef.2*coef.5  )

#cuadrantes$coef.1 + cuadrantes$coef.2*cuadrantes$coef.5

cuadrantes %>% group_by(allee.f,subexp.f) %>% 
  summarize_at(c("cociente","resta","angle"),.funs = list("median"=median,"iqr"=IQR),na.rm=T)

#filtrar casos con muy poco slope inicial y mucho final - cociente>5000
#cuadrantes=cuadrantes %>%  filter(cociente<5000,cociente>0)

# cuadrantes %>% ggplot(aes(x=log(cociente),fill=interaction(allee.f,subexp.f))) + geom_histogram(alpha=.3) + facet_wrap(~subexp.f*allee.f)
# cuadrantes %>% ggplot(aes(x=log(cociente),fill=allee)) + geom_histogram() 
# cuadrantes %>% ggplot(aes(x=resta,fill=allee.f)) + geom_histogram() + facet_grid(rows = "allee.f~subexp.f")
# 
# cuadrantes %>% ggplot(aes(x=log(cociente),fill=allee.f)) + geom_histogram() + facet_grid(rows = "allee.f~subexp.f",)
# 
# cuadrantes %>% ggplot(aes(y=log(cociente),x=log(1+coef.2),col=allee.f)) + geom_point(alpha=.2)+
#   facet_wrap(~subexp.f)+
#   labs(x="log(1+initial slope)",y="log(second slope/first slope)")
# 
# cuadrantes %>% ggplot(aes(y=coef.4,x=coef.2,col=allee.f)) + geom_point(alpha=.5)+
#   facet_wrap(~subexp.f)+
#   labs(x="log(first slope)",y="log(second slope)")
# 
cuadrantes$allee.f = factor(cuadrantes$allee.f)

# cuadrantes %>% mutate_at(c("coef.2","coef.3","cociente"),.funs = list("log"=function(x) log(1+x))) %>% 
# ggscatterhist(x="coef.2_log",y="cociente_log",
#   color="allee.f",alpha=0.6,size=3,
#   margin.plot="boxplot",
#   ggtheme=theme_bw(),
#   xlab="log(1+initial slope)",
#   ylab="log(1+final slope / initial slope)"
# ) 



##### 

colnames(cuadrantes)[c(10,12)]=c("initial.slope","slope.after.thresh")

cuadrantes=cuadrantes %>%  filter(initial.slope>0)

save(cuadrantes,file="simulados_allee_varia_parametros.RData")


##  plots varios con N
cuadrantes %>% ggplot(aes(x=log10(N),y=log10(I.active),col=allee.f))+geom_point() +facet_wrap(~allee.f)
cuadrantes %>% ggplot(aes(x=log10(N),y=log10(1+initial.slope),col=allee.f))+geom_point() +facet_wrap(~allee.f)



## -------------- PLOT CON CONVEX HULL -----------------
xlims=c(0,.5)
ylims= c(-.01,3)
plot_slopes =  # 
cuadrantes %>% mutate_at(vars(contains("slope"),"cociente"),.funs = list("log10"=function(x) log10(1+x))) %>%
 slice(sample(nrow(cuadrantes),1000)) %>% 
ggscatterhist(x="initial.slope_log10",y="cociente_log10",
              color="allee.f",alpha=.3,size=3,
              margin.plot="boxplot",
              ggtheme=theme_bw(),
              xlab="initial slope",
              ylab="slope after threshold",
              palette = c("#b33018","#14b74b"),
              ylim=ylims,
              xlim=xlims
)
plot_slopes$sp <-
plot_slopes$sp+
  geom_abline(slope=1,intercept = 0,linetype="dotted")# +  
 # geom_mark_hull(concavity = 5,expand = -.00000001,aes(fill=allee.f)) 


#plot_slopes$sp+
 # geom_abline(slope=1,intercept = 0,linetype="dotted") +  
#  geom_mark_hull(concavity = 5,expand = -.00000001,radius=.02,aes(fill=allee.f)) 



plot_slopes$sp$labels$colour=""
plot_slopes$sp$labels$fill=""
plot_slopes$yplot=plot_slopes$yplot + ylim(ylims)
plot_slopes$xplot=plot_slopes$xplot + ylim(xlims)
plot_slopes
  print_plot_slopes = print(plot_slopes)
ggsave(print_plot_slopes,filename = "fig1_slopes.pdf",height = 4,width = 5)

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


###################################################################
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


