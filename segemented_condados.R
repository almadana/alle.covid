#### Dara bases from The New York Times and from US Census
us.counties <- read.csv("us-counties.csv")
co.est2019.alldata <- read.csv("co-est2019-alldata.csv", header=TRUE)

nn<-cbind(as.numeric(paste(co.est2019.alldata[,4], 0,co.est2019.alldata[,5] , sep="")), co.est2019.alldata[,6:8]) # id or counties
colnames(nn)[1]<-"Id"

######### Function for ploting log10(infected) ~ log10(time) dynamics.  Infected =Io·(time)^mu following Maiers et al. 2020
condados<-function(nn, us.counties){
  library("segmented")
  out<-NULL
  #par(mfrow=c(5,5), mar=c(2,2,2,0.5 ))
  for (i in condPop[,"fips"]){
    ii<-which(us_counties[,"fips"]==i)
    if(length(ii)>15){
      x<-pull(us_counties[ii,5])
      X<-diff(x)
      id<-which(X==max(X))[1]
      if(length(14:id)>15 & min(x[14:id])>10 & max(x)>200) {    # 14 Days for having enough time for diseases dynamic (many counties have weeks with few cases)
        tt<-I(log10(1:length(x[14:id])))
        Infect<-log10(x[14:id])
        ff1<-segmented(obj = lm(Infect~ tt), seg.Z = ~tt, psi = 1)
        fit<-summary(ff1)
        p.seg<-coefficients(ff1)
        ff0<-lm(Infect~ tt)
        aic<-AIC(ff0,ff1)[,2]
        wi<-exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
  print("che")
        p.li<-coefficients(ff0)
 #   if(length(which(tt<fit$psi[2]))>10 & 
#       as.vector(p.seg)[2]<as.vector(p.seg)[3] &
#       as.vector(p.seg)[2]<1.2 &
#       length(which(tt>fit$psi[2]))>10 &
#       as.vector(p.seg)[3]>1.5
#       ){
        #plot(Infect~tt ,bty="l", pch=19, col="navy", cex=.8, main=nn[which(nn[,1]==i),3])        
        #plot(ff1, add=T, col="red")
#        }
        N.county<-nn[which(nn[,1]==i),4]                                          # Population of the country
        id.psi<-which(tt<fit$psi[2])                                              # The following lines estimate the number of active cases at the break point
        if(length(id.psi)>9) I.active<-Infect[max(id.psi)]-Infect[max(id.psi)-9]  # If the time series before the break is larger than 10 days retain the total infected 10 days before...
        if(length(id.psi)<=9)I.active<-Infect[max(id.psi)]                        # If the time series is shoreter take the total number of infected before the break
        out<-rbind(out,c(i,N.county, max(x),as.vector(p.seg)[-4],fit$psi[2],wi[2],fit$r.squared,10^I.active))
        ff0<-ff1<-NULL
      }
    }      
  }
  out<-cbind(out,10^(out[,4]+out[,5]*out[,7]))
  colnames(out)<-c("Id", "Population", "I.max", "intercept", "initial.slope", "slope.after.thresh", "time.threshold",
                   "Weight.evid","R.sqrt", "I.active", "Inf.total")
  out      
}      
save.image()
###########
condados(nn = condPop, us.counties = us_counties)->b

head(b)


b<-cbind(b,10^(b[,4]+b[,5]*b[,7]))
#b[,10]<-10^(b[,4]+b[,5]*b[,7])

#quitar los que no ajustgaron bien
b=  b[!(is.na(b[,8])),]

par(mfrow=c(2,3), mar=c(4,4,2,2))
hist(b[,5], col="navy", border = "gray", main="Initial Slope", xlab="Slope")

umbralBajo = b[,10]<10

hist(b[umbralBajo,6], col="red", border = "gray", main="Slope after  threshold", xlab="Slope",probability = T)
hist(b[!umbralBajo,6], col="navy", border = "gray", main="Slope after  threshold", xlab="Slope",add=T,probability = T)


hist(log10(b[,10]), col="navy", border = "gray", main="Active infected at breakpoint", xlab="Log10(Active Cases)")


#hist(log10(b[,10]), col="navy", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)

#hist((b[,10]/b[,2]*100), col="navy",
#     border = "gray", main="Infection Number ",
#     xlab="log10(Infected)", breaks=100, xlim=c(0,.5))

plot(b[,5]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Initial Slope")
abline(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
summary(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,5]~ log10(b[,2]),  col="gray")

plot(b[,6]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Slope after threshold")
#abline(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#summary(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,6]~ log10(b[,2]),  col="gray")

plot(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Log(Infection threshold)")
abline(lm(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2])), col="red", lwd=2)
summary(lm(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2])), col="red", lwd=2)
#text(4.5,3,"Inf.thres.~N^0.42\np:7.8e-12\nr-sqrt:0.20")
points(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2]),  col="gray")

####

plot(b[,5]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Initial Slope")
abline(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
summary(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,5]~ log10(b[,2]),  col="gray")

plot(b[,6]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Slope after threshold")
#abline(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#summary(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,6]~ log10(b[,2]),  col="gray")

plot(log10(b[,10])~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Log(Infection threshold)")
abline(lm(log10(b[,10])~ log10(b[,2])), col="red", lwd=2)
summary(lm(log10(b[,10])~ log10(b[,2])), col="red", lwd=2)
text(4.5,3,"Inf.thres.~N^0.42\np:7.8e-12\nr-sqrt:0.19")
points(log10(b[,10])~ log10(b[,2]),  col="gray")






#################
#paises










