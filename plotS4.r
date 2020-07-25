#makes plot if Supp Fig S4

par(mfrow=c(2,3), mar=c(4,6,2,2))
hist(b_p[,5], col="#4020ab", border = "gray", main="Initial Slope", xlab=expression(mu[1]))
mtext("a",side=3,at=-2.2,cex=1.6)
#text(1.25, 20, "b<1 potential growth rate\nlower than lineal increase")

hist((b_p[,6]+b_p[,5]), col="#4020ab", border = "gray", main="Slope after threshold", xlab=expression(mu[2]),10)
#text(6, 25, "b>1 potential growth rate\nhigher than lineal")
mtext("c",side=3,at=-6,cex=1.6)


hist(log10(b_p[,10]), col="#4020ab", border = "gray", main="Active infected at breakpoint ", xlab="Log10(Active Cases)", breaks=20)
mtext("e",side=3,at=-1.5,cex=1.6)
#hist((b_p[,10]/b_p[,2]*100), col="#4020ab",
#     border = "gray", main="Infection Number ",
#     xlab="log10(Infected)", breaks=100, xlim=c(0,.5))

plot(b_p[,5]~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab=expression(mu[1]))
abline(lm(b_p[,5]~ log10(b_p[,2])), col="red", lwd=2)
summary(lm(b_p[,5]~ log10(b_p[,2])), col="red", lwd=2)
#text(5.2,3,"b0~N^0.14\np:0.03\nr.sqr=0.04")
text(5.5,3,expression(mu[1] %~% Pop^0.14),cex=.8)
text(5.5,2.7,expression( p == .03),cex=.8)
text(5.5,2.4,expression(r^2 == 0.04),cex=.8)#,"\np:1.09e-12\nr.sq:0.3"))

points(b_p[,5]~ log10(b_p[,2]),  col="gray")
mtext("b",side=3,at=2.8,cex=1.6)

plot(b_p[,6]~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab=expression(mu[2]))

#abline(lm(b_p[,6]~ log10(b_p[,2])), col="red", lwd=2)
summary(lm(b_p[,6]~ log10(b_p[,2])), col="red", lwd=2)
#text(4.5,1.25,"b0~N^0.23\np:1.09e-07")
points(b_p[,6]~ log10(b_p[,2]),  col="gray")
mtext("d",side=3,at=2.8,cex=1.6)

plot(log10(b_p[,10])~ log10(b_p[,2]), bty="l", pch=19, col="#4020ab", xlab="Log10(Population)", ylab="Log(Active at threshold)")
#abline(lm(log10(b_p[b_p[,10]>10,10])~ log10(b_p[b_p[,10]>10,2])), col="red", lwd=2)
summary(lm(log10(b_p[b_p[,10]>10,10])~ log10(b_p[b_p[,10]>10,2])), col="red", lwd=2)
summary(lm(log10(b_p[b_p[,10]<10,10])~ log10(b_p[b_p[,10]<10,2])), col="red", lwd=2)
#text(4.5,3,"Active.thres.~N^0.19\np:7.8e-12\nr-sqrt:0.19")

points(log10(b_p[,10])~ log10(b_p[,2]),  col="gray")
mtext("f",side=3,at=2.8,cex=1.6)
