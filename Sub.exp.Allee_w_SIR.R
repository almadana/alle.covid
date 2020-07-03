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

par(mfrow=c(10,10), mar=c(2,2,0.5,0.5 ))
for(i in 1:100){
  dyn.sub.exp.Allee.SIR(I0=10, bet=0.8, gam=3 , p=0.8, t=100,  migration=1, size=0.2, I50=10, K=100000, Allee = F )->al
  imp<-cumsum(al[,5])  # cumulative number of imported cases
  al<-cumsum(al[,4])   # cumulative number of cases
  I.min<-which(al>9)
  al<-al[I.min]
  I.max<-which(al<0.8*1000)
  al<-al[I.max]
  print(length(al))
  al<-al[which(al>10 & al<2500)]
  imp<-imp[which(al>10 & al<2500)]
  if(length(al)>10){
    plot(log10(al)~log10(I(1:length(al))), bty="l", pch=19, col="darkgreen", main="")
  y<-log10(al)
  z<-log10(imp)
  x<-log10(1:length(y))
#  points(z~log10(I(1:length(al))), bty="l", pch=19, col="gold", main="", cex=.8)
  #ff0<-lm(y~x)
#  ff1<-segmented(lm(y~x), seg.Z = ~x, npsi = 1)
  #aic<-AIC(ff0,ffp, ff1)[,2]
  #wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
#  fit<-summary(ff1)
#  plot(ff1, add=T, col="red")
  }
}





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


