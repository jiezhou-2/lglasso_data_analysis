data=read.csv(file = "final.csv", header = T)[-c(248, 249),-1]
subject=data[,1]
a=data[,2]
age=rep(0,length(a))
for (i in 1:length(a)) {
age[i]=which(LETTERS==a[i])
}

rabundance=data[,3:ncol(data)]
lrdata=matrix(nrow = nrow(rabundance), ncol = ncol(rabundance)-1)
for (j in 1:ncol(lrdata)) {
  index=which(rabundance[,j]!=0)
  impute=min(rabundance[index,j])/10
  a1=ifelse(rabundance[,j]!=0,rabundance[,j],impute)
  a2=rabundance[,ncol(rabundance)]
  lrdata[,j]=log(a1/a2)
}


 lambda1=seq(0,1,length=1)
 lambda2=seq(1,5,length=1)
 lambda=c(lambda1,lambda2)
 #lambda=seq(0,1,length=5)
#ratio=seq(0,1,length=1)
ratio=20
 library(Matrix)
library(genlasso)
pool=modelpool(data = lrdata,age = age,group=group,lambda = lambda1,ratio = ratio)
pool=unique(pool)

data=as.matrix(data)
data1=lrdata[index1,]
data2=lrdata[index2,]
data3=lrdata[index3,]
data4=lrdata[index4,]
bicseq=rep(0,length(pool))
for (k in 1:length(pool)) {
  cat("_", k, "_")
  if (length(pool[[k]])==0) next()
theta1=pool[[k]][[1]]
theta2=pool[[k]][[2]]
theta3=pool[[k]][[3]]
theta4=pool[[k]][[4]]
prob1=0.2
prob2=0.2
prob3=0.2
prob4=0.2
P=0
bic1=gbic(data = data1,theta = theta1,prob = prob1,P=P)
bic2=gbic(data = data2,theta = theta2,prob = prob2,P=P)
bic3=gbic(data = data3,theta = theta3,prob = prob3,P=P)
bic4=gbic(data = data4,theta = theta4,prob = prob4,P=P)
bicseq[k]=bic1+bic2+bic3+bic4
}
k=which.min(bicseq)
bbinetwork=vector("list",length = 4)
family=colnames(data)[3:52]
networkmatrix=as.matrix(read.csv(file="network.csv")[,-1])
bbinetwork[[1]]=networkmatrix[,c(1:50)]
colnames(bbinetwork[[1]])=family
bbinetwork[[2]]=networkmatrix[,c(51:100)]
colnames(bbinetwork[[2]])=family
bbinetwork[[3]]=networkmatrix[,c(101:150)]
colnames(bbinetwork[[3]])=family
bbinetwork[[4]]=networkmatrix[,c(151:200)]
colnames(bbinetwork[[4]])=family

common12=as.matrix(bbinetwork[[1]]*bbinetwork[[2]])
common13=as.matrix(bbinetwork[[1]]*bbinetwork[[3]])
common14=as.matrix(bbinetwork[[1]]*bbinetwork[[4]])
common23=as.matrix(bbinetwork[[2]]*bbinetwork[[3]])
common24=as.matrix(bbinetwork[[2]]*bbinetwork[[4]])
common34=as.matrix(bbinetwork[[3]]*bbinetwork[[4]])
(sum(common12)-50)/2
(sum(common13)-50)/2
(sum(common14)-50)/2
(sum(common23)-50)/2
(sum(common24)-50)/2
(sum(common34)-50)/2
(sum(common)-50)/2
network1m2=ifelse(bbinetwork[[1]]==1 & bbinetwork[[2]]==0,1,0)
network1m3=ifelse(bbinetwork[[1]]==1 & bbinetwork[[3]]==0,1,0)
network1m4=ifelse(bbinetwork[[1]]==1 & bbinetwork[[4]]==0,1,0)
network2m3=ifelse(bbinetwork[[2]]==1 & bbinetwork[[3]]==0,1,0)
network2m4=ifelse(bbinetwork[[2]]==1 & bbinetwork[[4]]==0,1,0)
network3m4=ifelse(bbinetwork[[3]]==1 & bbinetwork[[4]]==0,1,0)

sum(network1m2)/2
sum(network1m3)/2
sum(network1m4)/2
sum(network2m3)/2
sum(network2m4)/2
sum(network3m4)/2


network2m1=ifelse(bbinetwork[[2]]==1 & bbinetwork[[1]]==0,1,0)
network3m1=ifelse(bbinetwork[[3]]==1 & bbinetwork[[1]]==0,1,0)
network4m1=ifelse(bbinetwork[[4]]==1 & bbinetwork[[1]]==0,1,0)
network3m2=ifelse(bbinetwork[[3]]==1 & bbinetwork[[2]]==0,1,0)
network4m2=ifelse(bbinetwork[[4]]==1 & bbinetwork[[2]]==0,1,0)
network4m3=ifelse(bbinetwork[[4]]==1 & bbinetwork[[3]]==0,1,0)

sum(network2m1)/2
sum(network3m1)/2
sum(network4m1)/2
sum(network3m1)/2
sum(network3m2)/2
sum(network4m3)/2



common=as.matrix(bbinetwork[[1]]*bbinetwork[[2]]*bbinetwork[[3]]*bbinetwork[[4]])
colnames(common)=family
library(igraph)
library(dplyr)
library(ggplot2)
for (k in 1:4) {
  N=common
  #colnames(M)=paste("M",1:50,sep = "")
  g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "directed")
  #V(g2)$comm=membership(optimal.community(g2))
  V(g2)$degree = degree(g2)
  #V(g2)$closeness = centralization.closeness(g2)$res
  #V(g2)$betweenness = centralization.betweenness(g2)$res
  #V(g2)$eigen = centralization.evcent(g2)$vector
  #V(g2)$name=paste("M",1:50, sep = "")
  V(g2)$name=family
  node_list=get.data.frame(g2,what="vertices")
  edge_list=get.data.frame(g2,what="edges")
  group=rep(1,nrow(edge_list))
  edge_list=cbind(edge_list,group)
  all_nodes <- sort(node_list$name)
  plot_data <- edge_list %>% mutate(
    to = factor(to, levels = all_nodes),
    from = factor(from, levels = all_nodes))
  png(filename = paste("networkcommon", ".png", sep = ""))
  ggplot(plot_data, aes(x = from, y = to, fill = group,xlab="Microbe",ylab="Microbes"))  +
    geom_raster() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "none")
  dev.off()
}




