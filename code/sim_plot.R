load("./data/homoc2/uu2=0.3.Rd")
results01=simures
load("./data/homoc2/uu2=1.3.Rd")
results02=simures
load("./data/homoc2/uu2=2.3.Rd")
results03=simures
load("./data/homoc2/uu2=3.3.Rd")
results04=simures



par(mfrow=c(2,2),mar=c(4,4,2,2),oma=c(0,0,2,0))
FPR=results01[[1]][[1]][,2]
TPR=results01[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results01[[1]][[2]][,2],results01[[1]][[2]][,1],type="l",lty=2)
lines(results01[[1]][[3]][,2],results01[[1]][[3]][,1],type="l",lty=3)
legend(0.5,0.4,legend=c("(2.3,0.3)"))


FPR=results02[[1]][[1]][,2]
TPR=results02[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results02[[1]][[2]][,2],results02[[1]][[2]][,1],type="l",lty=2)
lines(results02[[1]][[3]][,2],results02[[1]][[3]][,1],type="l",lty=3)
legend(0.5,0.4,legend=c("(2.3,1.3)"))





FPR=results03[[1]][[1]][,2]
TPR=results03[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results03[[1]][[2]][,2],results03[[1]][[2]][,1],type="l",lty=2)
lines(results03[[1]][[3]][,2],results03[[1]][[3]][,1],type="l",lty=3)
legend(0.5,0.4,legend=c("(2.3,2.3)"))

FPR=results04[[1]][[1]][,2]
TPR=results04[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results04[[1]][[2]][,2],results04[[1]][[2]][,1],type="l",lty=2)
lines(results04[[1]][[3]][,2],results04[[1]][[3]][,1],type="l",lty=3)
legend(0.5,0.4,legend=c("(2.3,3.3)"))


title("Homogeneous Subjects with Heterogeneous Microbial Communities", outer=TRUE, cex=1.5)

