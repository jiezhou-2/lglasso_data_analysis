##compare the correlation of distance between precision network and phylogenetic network
distance2=distances(treegraph)
index=c()
for (i in 1:length(family)) {
  index=c(index,which(rownames(distance2)==family[i]))
}

subdistance2=distance2[index,index]
###
rrr=c()
D=vector("list",length = length(pathlist))
for (i in 1:length(pathlist)) {
  N=pathlist[[i]]
  diag(N)=0
  g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "undirected",diag = F)  
  distance1=distances(g2)
  index=which(distance1==Inf,arr.ind = T)
  distance1[index]=100
  d1=distance1[lower.tri(distance1)]
  d2=subdistance2[lower.tri(subdistance2)]
  d1_NM=d1[which(d1!=100)]
  d2_NM=d2[which(d1!=100)]
  rrr=c(rrr,cor(d1_NM,d2_NM))
  d_NM=cbind(d1_NM,d2_NM)
  D[[i]]=d_NM
  print(length(d1_NM))
}
plot(rrr,type="l",ylim=c(0,0.5))
lines(u[30:63]/3321)

