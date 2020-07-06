library(segmented)


casos.world.264 <- read.delim("./casos.world.264.tab")

par(mfrow=c(5,15), mar=c(1,1,2,0.5 ))
paises<-function(casos.world.264){
  library("segmented")
  out<-NULL
  for (i in 1:nrow(casos.world.264)){
    ii<-which(casos.world.264[i,18:154]>=1)
    if(length(ii)>15){
      x<-as.numeric(casos.world.264[i,(ii[1]+17):154])
      X<-diff(x)
      id<-which(X==max(X))[1]
      #id<-which(x==max(x))[1]
      if(length(14:id)>15 & min(x[14:id])>20 & max(x)>100) {
        tt<-I(log10(1:length(x[14:id])))
        Infect<-log10(x[14:id])
        ff0<-lm(Infect~ tt)
        ff1<-segmented(obj = lm(Infect~ tt), seg.Z = ~tt, npsi = 1)
        aic<-AIC(ff0,ff1)[,2]
        wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
        N.pais<-casos.world.264[i,13]
        p.li<-coefficients(ff0)
        p.seg<-coefficients(ff1)
        fit<-summary(ff1)
       # if(length(which(tt<fit$psi[2]))>9 & 
      #     as.vector(p.seg)[2]<as.vector(p.seg)[3] &
      #     as.vector(p.seg)[2]<1.1 &
      #     length(which(tt>fit$psi[2]))>10 &
      #     as.vector(p.seg)[3]>1.1
      #  ){
        plot(Infect~tt ,bty="l", pch=19, col="navy", cex=.8, main=casos.world.264[i,9]) 
        plot(ff1, add=T, col="red")
      #x  }
  #      bbreak<-which(x>fit$psi[2])[1]
        
        out<-rbind(out,c(i,N.pais, max(x),as.vector(p.seg)[-4],fit$psi[2],wi[2],fit$r.squared))
      }
    }      
  }
  colnames(out)<-c("Id", "Population", "I.max", "intercept", "initial.slope", "slope.after.thresh", "time.threshold","Weight.evid","R.sqrt")
  out      
}      
save.image()
###########
paises(casos.world.264)->b
head(b)
dim(b)
b<-cbind(b,10^(b[,4]+b[,5]*b[,7]))
#b[,10]<-10^(b[,4]+b[,5]*b[,7])
par(mfrow=c(2,3), mar=c(4,4,2,2))
hist(b[,5], col="navy", border = "gray", main="Initial Slope", xlab="Slope")
text(1.25, 20, "b<1 potential growth rate\nlower than lineal increase")

hist((b[,6]+b[,5]), col="navy", border = "gray", main="Slope After Threshold", xlab="Slope")
text(6, 25, "b>1 potential growth rate\nhigher than lineal")

hist(log10(b[,10]), col="navy", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)

#hist((b[,10]/b[,2]*100), col="navy",
#     border = "gray", main="Infection Number ",
#     xlab="log10(Infected)", breaks=100, xlim=c(0,.5))

plot(b[,5]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log10(Population)", ylab="Initial Slope")
abline(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
summary(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,5]~ log10(b[,2]),  col="gray")

plot(b[,6]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log10(Population)", ylab="Slope after threshold")
#abline(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
summary(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,6]~ log10(b[,2]),  col="gray")

plot(log10(b[,10])~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log10(Population)", ylab="Log(Infection threshold)")
abline(lm(log10(b[,10])~ log10(b[,2])), col="red", lwd=2)
summary(lm(log10(b[,10])~ log10(b[,2])), col="red", lwd=2)
text(4.5,3,"Inf.thres.~N^0.42\np:7.8e-12\nr-sqrt:0.19")
points(log10(b[,10])~ log10(b[,2]),  col="gray")

####

#################
#paises

colnames(b)
hist((b[,8]), col="navy", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)
hist((b[,9]), col="navy", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)






