
##beta-diversity
library(dada2)
library(DECIPHER)
library(phangorn)
library(phytools)
library(igraph)
library(ggplot2)


##sequence bootstrap
library(boot)
fc <- function(d, i){
  dd <- d[i,]
  return(cor(dd[,1], dd[,2]))
}
aa=c()
bb=c()
i=33
for (i in 1:40) {
bootcor=boot(D[[i]],fc,R=5000)
aa=c(aa,mean(bootcor$t))
bb=c(bb,sd(bootcor$t))
}

hist(bootcor$t,breaks = 200,xlab = "Correlation for bootstrapped distances",main = "")

