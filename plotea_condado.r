plot_fit_condado <- function(i) {
  require(segmented)
  ii<-which(us.counties[,4]==i)
  if(length(ii)>15){
    x<-us.counties[ii,5]
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
      
      p.li<-coefficients(ff0)
      #   if(length(which(tt<fit$psi[2]))>10 & 
      #       as.vector(p.seg)[2]<as.vector(p.seg)[3] &
      #       as.vector(p.seg)[2]<1.2 &
      #       length(which(tt>fit$psi[2]))>10 &
      #       as.vector(p.seg)[3]>1.5
      #       ){
      plot(Infect~tt ,bty="l", pch=19, col="navy", cex=.8, main=nn[which(nn[,1]==i),3])        
      plot(ff1, add=T, col="red")
      #        }
      N.county<-nn[which(nn[,1]==i),4]                                          # Population of the country
      id.psi<-which(tt<fit$psi[2])                                              # The following lines estimate the number of active cases at the break point
      if(length(id.psi)>9) I.active<-Infect[max(id.psi)]-Infect[max(id.psi)-9]  # If the time series before the break is larger than 10 days retain the total infected 10 days before...
      if(length(id.psi)<=9)I.active<-Infect[max(id.psi)]                        # If the time series is shoreter take the total number of infected before the break
      
      return(ff1)
    }
  }
}

