casos.world.264 <- read.delim("casos.world.264.tab")


paises<-function(casos.world.264,nCols){
  library("segmented")
  out<-NULL
  par(mfrow=c(10,10), mar=c(2,2,3,0.5 ))
  for (i in 1:nrow(casos.world.264)){
#    ii<-which(casos.world.264[i,18:154]>=1) #datos viejos
    ii<-which(casos.world.264[i,18:nCols]>=1) #al 10 de julio
    if(length(ii)>15){
      x<-as.numeric(casos.world.264[i,(ii[1]+17):nCols])
      X<-diff(x)
      id<-which(X==max(X))[1]
      #id<-which(x==max(x))[1]
      if(length(14:id)>15 & min(x[14:id])>9 & max(x)>100) {
        tt<-I(log10(1:length(x[14:id])))
        Infect<-log10(x[14:id])
        #plot(Infect~tt ,bty="l", pch=19, col="#4020ab", cex=.8, main=casos.world.264[i,9])        
        ff0<-lm(Infect~ tt)
        ff1<-segmented(obj = lm(Infect~ tt), seg.Z = ~tt, npsi = 1)
        aic<-AIC(ff0,ff1)[,2]
        wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
        N.pais<-casos.world.264[i,"Population"]
        p.li<-coefficients(ff0)
        p.seg<-coefficients(ff1)
        fit<-summary(ff1)
        #plot(ff1, add=T, col="red")
        bbreak<-which(x>fit$psi[2])[1]
        id.psi<-which(tt<fit$psi[2])                                              # The following lines estimate the number of active cases at the break point
        if(length(id.psi)>9) I.active<-Infect[max(id.psi)]-Infect[max(id.psi)-9]  # If the time series before the break is larger than 10 days retain the total infected 10 days before...
        if(length(id.psi)<=9)I.active<-Infect[max(id.psi)] 
        
        
        out<-rbind(out,c(i,N.pais, max(x),as.vector(p.seg)[-4],fit$psi[2],wi[2],fit$r.squared,10^I.active))
      }
    }      
  }
  colnames(out)<-c("Id", "Population", "I.max", "intercept", "initial.slope", "slope.after.thresh", "time.threshold","Weight.evid","R.sqrt","I.active")
  out      
}      
save.image()
###########
#paises(casos.world.264)->b
paises(cases.world,188)->b_p
head(b_p)


b_p<-cbind(b_p,10^(b_p[,4]+b_p[,5]*b_p[,7]))
colnames(b_p)[11]="I.total"

#quitar Eritrea

b_p=b_p[-40,]
#
#b_p[,10]<-10^(b_p[,4]+b_p[,5]*b_p[,7])
par(mfrow=c(2,3), mar=c(4,4,2,2))
hist(b_p[,5], col="#4020ab", border = "gray", main="Initial Slope", xlab="Slope")
#text(1.25, 20, "b<1 potential growth rate\nlower than lineal increase")

hist((b_p[,6]+b_p[,5]), col="#4020ab", border = "gray", main="Slope after threshold", xlab="Slope",10)
#text(6, 25, "b>1 potential growth rate\nhigher than lineal")

hist(log10(b_p[,10]), col="#4020ab", border = "gray", main="Active infected at breakpoint ", xlab="Log10(Active Cases)", breaks=20)

#hist((b_p[,10]/b_p[,2]*100), col="#4020ab",
#     border = "gray", main="Infection Number ",
#     xlab="log10(Infected)", breaks=100, xlim=c(0,.5))

plot(b_p[,5]~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab="Initial Slope")
abline(lm(b_p[,5]~ log10(b_p[,2])), col="red", lwd=2)
summary(lm(b_p[,5]~ log10(b_p[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b_p[,5]~ log10(b_p[,2]),  col="gray")

plot(b_p[,6]~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab="Slope after threshold")

#abline(lm(b_p[,6]~ log10(b_p[,2])), col="red", lwd=2)
summary(lm(b_p[,6]~ log10(b_p[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b_p[,6]~ log10(b_p[,2]),  col="gray")

plot(log10(b_p[,10])~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab="Log(Infection threshold)")
abline(lm(log10(b_p[,10])~ log10(b_p[,2])), col="red", lwd=2)
summary(lm(log10(b_p[,10])~ log10(b_p[,2])), col="red", lwd=2)
#text(4.5,3,"Inf.thres.~N^0.42\np:7.8e-12\nr-sqrt:0.19")

points(log10(b_p[,10])~ log10(b_p[,2]),  col="gray")

####

#################
#paises

colnames(b_p)
hist((b_p[,8]), col="#4020ab", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)
hist((b_p[,9]), col="#4020ab", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)





