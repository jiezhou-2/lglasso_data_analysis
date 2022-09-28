

##Define the testing statistic
theta=function(v1,v2){
  r1=cor(v1,v2)
  v12=var(v1)
  v22=var(v2)
  v12v22=mean((v1-mean(v1))^2*(v2-mean(v2))^2)
  vv=v12v22/(v12*v22)
  ss=sqrt(length(v1))*r1/sqrt(vv)
  return(ss)
}


##sequence simulation
d2=subdistance2[lower.tri(subdistance2)]
L=matrix(nrow = 5001,ncol = 40)
for (i in 1:length(pathlist)) {
  print(i)
  g=graph_from_adjacency_matrix(adjmatrix = t(pathlist[[i]]), mode = "undirected",diag = F)
  distance1=distances(g)
  d1=distance1[lower.tri(distance1)]
  d1_NM=d1[which(d1!=Inf)]
  d2_NM=d2[which(d1!=Inf)]
  rr=theta(d1_NM,d2_NM)
  for (j in 1:5000) {
  g2=permute.vertices(g,sample(vcount(g)))
  distance1=distances(g2)
  d1=distance1[lower.tri(distance1)]
  d1_NM=d1[which(d1!=Inf)]
  d2_NM=d2[which(d1!=Inf)]
  rr=c(rr,theta(d1_NM,d2_NM))
  }
  L[,i]=rr
}
##plot
boxplot(L,ylab="Studentized Correlation", xlab="Model Index",ylim=c(-17,17))
lines(L[1,])

boxplot(L[-1,],ylab="Studentized Correlation", xlab="Model Index")


##point permutation
N=pathlist[[43]]
diag(N)=0
g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "undirected",diag = F)  
distance1=distances(g2)
index=which(distance1==Inf,arr.ind = T)
distance1[index]=100
d1=distance1[lower.tri(distance1)]
d2=subdistance2[lower.tri(subdistance2)]
d1_NM=d1[which(d1!=100)]
d2_NM=d2[which(d1!=100)]
D=cbind(d1_NM,d2_NM)
cor(d1_NM,d2_NM)
r0=theta(v1=d1_NM,d2_NM)
per=5000
pcor=c()
pcor1=c()
for (i in 1:per) {
  print(i)
  g22=permute.vertices(g2,sample(vcount(g2)))
  distance1=distances(g22)
  index=which(distance1==Inf,arr.ind = T)
  distance1[index]=100
  d1=distance1[lower.tri(distance1)]
  d2=subdistance2[lower.tri(subdistance2)]
  d1_NM=d1[which(d1!=100)]
  d2_NM=d2[which(d1!=100)]
  pcor=c(pcor,cor(d1_NM,d2_NM))
  pcor1=c(pcor1,theta(d1_NM,d2_NM))
  
}
hist(pcor,breaks = 20,xlab = "Correlation",main = "", xlim = c(-0.3,0.3))
abline(v= 0.2284483,lty=3,col=2)
text(0.25,0.1,"   observed correlation")
mean(pcor)
sd(pcor)
hist(pcor1,breaks = 200,xlab = "Permuted studentized correlation",main = "",xlim = c(-15,15))
abline(v= r0,lty=3,col=2)
text(13,0.1,"   observed")
mean(pcor1)
sd(pcor1)


##point bootstrap
library(boot)
fc <- function(d, i){
  dd <- d[i,]
  return(cor(dd[,1], dd[,2]))
}

  bootcor=boot(D,fc,R=5000)
  ci=quantile(bootcor$t,probs = c(0.025,0.975))
  hist(bootcor$t,breaks = 100,xlab = "Correlation",main="")
  abline(v=  0.2440453,lty=3,col=2)
  text( 0.26,-1,"   observed correlation")

  realcorrboot=data.frame(x=bootcor$t)
  theme_set(theme_classic()+
              theme(legend.position = "top"))
  e=ggplot(realcorrboot,aes(x=x))
  e+geom_histogram(aes(y=stat(density)),binwidth=0.001, colour="black",fill="white")+
    #geom_density(alpha=0.2,fill="#E7B800")+
    geom_vline(aes(xintercept=mean(x)),linetype="dashed",size=0.6,colour="red")+
    labs(x="Correlation",y="Density")
  
  
  mean(bootcor$t)
sd(bootcor$t)




#index2=which(distance1==Inf,arr.ind = T)
#distance1[index2]=100
# d1=distance1[lower.tri(distance1)]
# d2=subdistance2[lower.tri(subdistance2)]
# d1_NM=d1[which(d1!=Inf)]
# d2_NM=d2[which(d1!=Inf)]
# rr=c(rr,cor(d1_NM,d2_NM))
# for (i in 1:per) {
#   perfamily=family[sample()]
#   perd1=d1[sample(length(82))]
#   pcor=c(pcor,cor(perd1,d2))
# }
# hist(pcor,breaks=200, xlab = "Correlation for permuted data",
#      xlim = c(-0.07,0.15),main ="")
# 
# abline(v=cor(d1,d2),lty=3,col=2)
# text(0.12,0.1,"   observed correlation")


###compare the intra-community distance and inter-community distance
g=function(network,community){
ad=as_adjacency_matrix(network)
nodes=rownames(ad)
p=length(nodes)
m=length(community)
index=c()
for (i in 1:length(community)) {
  index=c(index,which(nodes==community[i]))
}
supplement=setdiff(nodes,community)
intradismat=distances(network,v=community,to=community)
interdismat=distances(network,v=community,to=supplement)

intradis=intradismat[lower.tri(intradismat)]
interdis=c(interdismat)
d=matrix(nrow = 2,ncol = 3)
d[1,1]=mean(intradis)
d[1,2]=mean(interdis)
d[1,3]=d[1,1]-d[1,2]
d[2,1]=sd(intradis)
d[2,2]=sd(interdis)
d[2,3]=sqrt(d[2,1]^2/(m*(m-1)/2)+d[2,2]^2/(m*p))
return(d)
  }

test1=g(network = treegraph,community = commu1)
test2=g(network = treegraph,community = commu2)
test3=g(network = treegraph,community = commu3)
test4=g(network = treegraph,community = commu4)
test5=g(network = treegraph,community = commu5)
test6=g(network = treegraph,community = commu6)
test7=g(network = treegraph,community = commu7)
test8=g(network = treegraph,community = commu8)
# 
# graphlanBigTree=data.frame(matrix(nrow = nrow(taxa),ncol = 1))
# for (i in 1:nrow(taxa)){
#   index=length(which(!is.na(taxa[i,])))
#   a=taxa[i,1]
#   if (index>1){
#     for (j in 2:index) {
#       a=paste(a,taxa[i,j],sep = ".")
#     }
#   }
#   
#   graphlanBigTree[i,]=a
# }
# write.csv(graphlanBigTree,file="graphlanBigTree.csv")

