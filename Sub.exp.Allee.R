



###################################################################
##### Dynamic

pred.dyn.sub.exp.Allee<-function(I0, r, p, K, t, Allee=T, migration, size=0.5, I50=50, f0=0.3){
  out<-I0
  Inf.active<-c(rep(0,6),I0)                                                      # This will record the number of active cases. It is used for the evaluation of the Allee Effect
  for(i in 1:t){
    It1<-out[length(out)]
    if(It1<=0){
      It2<-0
    } else{
    active.cases<-sum(Inf.active)
    if(Allee==T)D.I<-(r*It1^p)*((1/(1-f0))*(active.cases/(active.cases+I50))-f0)*(1-It1/K)
    if(Allee==F)D.I<-(r*It1^p)*(1-It1/K)
    imported<-rnbinom(n = 1, size = size, mu = migration)
    if(D.I>0)D.I<-rnbinom(n = 1, size = size, mu = D.I)+imported
    if(D.I<=0)D.I<-D.I+imported
    It2<-It1+D.I
    Inf.active<-c(Inf.active[-1],D.I)
    }
    if(It2<=0)It2<-0
#  print(c(It1, It2, D.I, imported))
    out<-c(out,It2)
  }
#  print(out)
  out
}



###################################################################
##### plot with/wihout Allee and subexponential

plot.all.su.exp<-function(I0 , r , p , K,t, migration, size, I50=50, f0=0.3, it){
require(segmented)
  par(mfrow=c(2,2), mar=c(4,4,4,4))
  for(pp in c(1,p)){
    for(m in c(0,migration)){
      for(allee in c(F,T)){
    out.t<-NULL
    for(i in 1:it){
    pred.dyn.sub.exp.Allee(I0 = I0, r = r, p = pp, K=K,t = t,
                  Allee = allee, migration=m, size=size, I50=50, f0=0.3)->al
#      al<-cumsum(al)
      out.t<-rbind(out.t, al)
      }
out<- apply(out.t, 2, quantile, c(0.025,0.5,0.975))
#    ii<-which(diff(al)==max(diff(al)))
      ii<-t
      Total.Infected<-out[2,]
      if(allee ==F & pp==1 & m ==0)         title<-"Without Allee & Without Sub.exp without migra"
      if(allee ==F & pp<1 & m ==0)          title<-"Without Allee & With Sub.exp without migra"
      if(allee ==T & pp<1 & m ==0)          title<-"With Allee & With Sub.exp without migra"
      if(allee ==T & pp==1 & m ==0)         title<-"With Allee & Without Sub.exp without migra"
      if(allee ==F & pp==1 & m ==migration) title<-"Without Allee & Without Sub.exp  with migra"
      if(allee ==F & pp<1 & m ==migration)  title<-"Without Allee & With Sub.exp  with migra"
      if(allee ==T & pp<1 & m ==migration)  title<-"With Allee & With Sub.exp   with migra"
      if(allee ==T & pp==1 & m ==migration) title<-"With Allee & Without Sub.exp   with migra"
#return(Total.Infected)      
    if(allee==F){  plot(log10(Total.Infected)~log10(I(1:length(Total.Infected))), bty="l", pch=19, col="navy", main=title, ylim=c(0,max(log10(out), na.rm=T)))
      points(log10(out[1,])~log10(I(1:length(Total.Infected))), col="navy", lty=2, type="l")
      points(log10(out[3,])~log10(I(1:length(Total.Infected))), col="navy", lty=2, type="l")
      }
      if(allee==T){
        points(log10(Total.Infected)~log10(I(1:length(Total.Infected))), bty="l", pch=19, col="darkred", main=title, ylim=c(0,max(log10(out), na.rm=T)))
        points(log10(out[1,])~log10(I(1:length(Total.Infected))), col="darkred", lty=2, type="l")
        points(log10(out[3,])~log10(I(1:length(Total.Infected))), col="darkred", lty=2, type="l")
        
          abline(h=log10(I50/((1/f0)-1)), col="blue", lty=2)      # This is the epidemic threshold because Allee effect
print(c(log10(I50/((1/f0)-1)),(I50/((1/f0)-1))))        
          y<-log10(Total.Infected)
          x<-log10(1:length(y))
        #ff0<-lm(y~x)
        #ffp<-lm(y~x+I(x^ 2))
          ff1<-segmented(lm(y~x), seg.Z = ~x, npsi = 1)
        #aic<-AIC(ff0,ffp, ff1)[,2]
        #wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
        #text(quantile(x,0.2),quantile(y,0.8),paste("wi",round(wi[3],3), paste=""), pos=4)
          fit<-summary(ff1)
          plot(ff1, add=T, col="red")
      #text(quantile(x,0.2),quantile(y,0.7),paste("r-sqrt:",round(fit$r.squared,3), paste=""), pos=4)
        }
      }
    }
  }
}


dd<-plot.all.su.exp(I0 = 10, r = 0.25, p = 0.75, K = 20000, t = 50, migration=1, size=.5, I50 = 50, f0 = 0.3, it=200)

plot.all.su.exp(I0 = 80, r = 0.15, p = 0.9, thr = 50, K = 1000, t = 50,  migration = 10, size = 1, Q=0.85)


l<-pred.dyn.sub.exp.Allee(I0 = 10, r = 0.25, p = 0.85, K = 450, t = 100, Allee = T, migration = F, size = 10, I50 = 50, f0 = 0.3)
