##compare the correlation of distance between precision network and phylogenetic network
distance2=distances(treegraph)
index=c()
for (i in 1:length(family)) {
  index=c(index,which(rownames(distance2)==family[i]))
}

subdistance2=distance2[index,index]
###
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso_data_analysis/data/network1.rd")
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso_data_analysis/data/network2.rd")

load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso_data_analysis/data/network3.rd")

rrr=matrix(nrow=3,ncol = length(network1))
yyy=matrix(nrow=3,ncol = length(network1))
pathlist=network1
D=vector("list",length = length(pathlist))
for (i in 1:length(pathlist)) {
  N1=ifelse(network1[[i]]==0,0,1)
  N2=ifelse(network2[[i]]==0,0,1)
  N3=ifelse(network3[[i]]==0,0,1)
  diag(N1)=0
  diag(N2)=0
  diag(N3)=0
  g1=graph_from_adjacency_matrix(adjmatrix = t(N1), mode = "undirected",diag = F)
  distance1=distances(g1)

  g2=graph_from_adjacency_matrix(adjmatrix = t(N2), mode = "undirected",diag = F)
  distance2=distances(g2)

  g3=graph_from_adjacency_matrix(adjmatrix = t(N3), mode = "undirected",diag = F)
  distance3=distances(g3)


  yyy[1,i]=sum(N1)/2
  yyy[2,i]=sum(N2)/2
  yyy[3,i]=sum(N3)/2

  index1=which(distance1==Inf,arr.ind = T)
  index2=which(distance2==Inf,arr.ind = T)
  index3=which(distance3==Inf,arr.ind = T)


  distance1[index1]=100
  distance2[index2]=100
  distance3[index3]=100


  d1=distance1[lower.tri(distance1)]
  d2=distance2[lower.tri(distance2)]
  d3=distance3[lower.tri(distance3)]
  d4=subdistance2[lower.tri(subdistance2)]


  d1_NM=d1[which(d1!=100)]
  d4_NM=d4[which(d1!=100)]
  rrr[1,i]=cor(d1_NM,d4_NM)

  d2_NM=d2[which(d2!=100)]
  d4_NM=d4[which(d2!=100)]
  rrr[2,i]=cor(d2_NM,d4_NM)

  d3_NM=d3[which(d3!=100)]
  d4_NM=d4[which(d3!=100)]
  rrr[3,i]=cor(d3_NM,d4_NM)


  print(length(d1_NM))
  print(length(d2_NM))
  print(length(d3_NM))
}
lgl=cbind(x=yyy[1,],y=rrr[1,])
lgl=lgl[order(lgl[,1]),]
gl=cbind(x=yyy[2,],y=rrr[2,])
gl=gl[order(gl[,1]),]
nh=cbind(x=yyy[3,],y=rrr[3,])
nh=nh[order(nh[,1]),]
plot(lgl[,1],lgl[,2],type="l",lty=1, ylim=c(0,0.55),xlim = c(20,1400),
     xlab = "Number  of edges ", ylab = "Correlation")
lines(gl[,1],gl[,2],lty=2)
lines(nh[,1],nh[,2],lty=3)
legend(900,0.55,legend = c("LGLASSO","GLASSO","NH"),lty=c(1,2,3))


