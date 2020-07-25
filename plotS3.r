#makes plot if Supp Fig S3

par(mfrow=c(2,3), mar=c(4,6,2,2))
hist(b[,5], col="navy", border = "gray", main="Initial Slope", xlab=expression(mu[1]))
mtext("a",side=3,at=-.8,cex=1.6)


hist(b[,6], col="navy", border = "gray", main="Slope after threshold", xlab=expression(mu[2]))
mtext("c",side=3,at=-6,cex=1.6)


hist(log10(b[,10]), col="navy", border = "gray", main="Active infected at breakpoint", xlab="Log10(Active Cases)")
mtext("e",side=3,at=-1.4,cex=1.6)


#hist(log10(b[,10]), col="navy", border = "gray", main="Infection Threshold ", xlab="log10(Infected)", breaks=20)

#hist((b[,10]/b[,2]*100), col="navy",
#     border = "gray", main="Infection Number ",
#     xlab="log10(Infected)", breaks=100, xlim=c(0,.5))

plot(b[,5]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab=expression(mu[1]))
abline(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
summary(lm(b[,5]~ log10(b[,2])), col="red", lwd=2)
#text(4.8,1.5,paste0(expression(mu[1] ~ N^0.33),"\np:1.09e-12\nr.sq:0.3"))
text(4.8,1.5,expression(mu[1] %~% Pop^0.33),cex = 0.8)
text(4.8,1.3,expression( p<10^12),cex = 0.8)
text(4.8,1.1,expression(r^2 == 0.3),cex = 0.8)#,"\np:1.09e-12\nr.sq:0.3"))

points(b[,5]~ log10(b[,2]),  col="gray")
mtext("b",side=3,at=3,cex=1.6)



plot(b[,6]~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab=expression(mu[2]))
#abline(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#summary(lm(b[,6]~ log10(b[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b[,6]~ log10(b[,2]),  col="gray")
mtext("d",side=3,at=3,cex=1.6)


#plot(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Log(Infection threshold)")
plot(log10(b[,10])~ log10(b[,2]), bty="l", pch=19, col="navy", xlab="Log(Population)", ylab="Log(Infection threshold)")
abline(lm(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2])), col="red", lwd=2)
summary(lm(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2])), col="red", lwd=2)
#text(5,3,"Inf.thres.~N^0.61\np:7.8e-10\nr.sq:0.34")
text(5,3,expression("Inf. thres." %~% Pop^0.61),cex = 0.8)
text(5,2.7,expression( p<10^9),cex = 0.8)
text(5,2.4,expression(r^2 == 0.34),cex = 0.8)#,"\np:1.09e-12\nr.sq:0.3"))

#points(log10(b[-which(log10(b[,10])<1),10])~ log10(b[-which(log10(b[,10])<1),2]),  col="gray")
points(log10(b[,10])~ log10(b[,2]),  col="gray")
mtext("f",side=3,at=3,cex=1.6)
