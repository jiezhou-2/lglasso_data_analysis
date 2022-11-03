data1=read.csv(file = "./data/final.csv", header = T)[,-1]
subject=data[,1]
age=data[,2]
rabundance=data[,3:ncol(data)]
lrdata=matrix(nrow = nrow(rabundance), ncol = ncol(rabundance)-1)
for (j in 1:ncol(lrdata)) {
  index=which(rabundance[,j]!=0)
impute=min(rabundance[index,j])/10
    a1=ifelse(rabundance[,j]!=0,rabundance[,j],impute)
  a2=rabundance[,ncol(rabundance)]
  lrdata[,j]=log(a1/a2)
}
#network=lnetwork(data = lrdata, subject=subject,age=age )
network=read.csv(file = "fullnetwork.csv")[,-1]
binetwork=ifelse(abs(network)<10^(-4),0,1)
bbinetwork=vector("list", 19)
for (i in 1:19) {
  Mi=binetwork[1:50,(1+49*(i-1)):(49*i)]
  M=matrix(0,nrow = nrow(Mi), ncol = ncol(Mi)+1)
  for (j in 1:nrow(M)) {
    M[j,-j]=Mi[j,]
  }
  #entry-wise multiplication
  M=M*t(M)
  bbinetwork[[i]]=M
}
library(igraph)
library(dplyr)
library(ggplot2)
for (k in 1:19) {
N=bbinetwork[[k]]
#colnames(M)=paste("M",1:50,sep = "")
g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "directed")
#V(g2)$comm=membership(optimal.community(g2))
V(g2)$degree = degree(g2)
#V(g2)$closeness = centralization.closeness(g2)$res
#V(g2)$betweenness = centralization.betweenness(g2)$res
#V(g2)$eigen = centralization.evcent(g2)$vector
V(g2)$name=paste("M",1:50, sep = "")
node_list=get.data.frame(g2,what="vertices")
edge_list=get.data.frame(g2,what="edges")
group=rep(1,nrow(edge_list))
edge_list=cbind(edge_list,group)
all_nodes <- sort(node_list$name)
plot_data <- edge_list %>% mutate(
  to = factor(to, levels = all_nodes),
  from = factor(from, levels = all_nodes))
png(filename = paste("graph",k, ".png", sep = ""))
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
#plot.igraph(g2,layout=layout_with_lgl, vertex.color="green", vertex.size=5)
#write.csv(network,file = "fullnetwork.csv")



