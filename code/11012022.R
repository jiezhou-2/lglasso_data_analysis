#estnetwork1=(estnetwork1+t(estnetwork1))
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso_data_analysis/data/estnetwork1.rd")
estnetwork1=10*(estnetwork1+t(estnetwork1))
diag(estnetwork1)=0
family=colnames(data)[-c(1,2)]
rownames(estnetwork1)=family
colnames(estnetwork1)=family
family1=family
family2=family
l=length(family)
for (i in 1:l) {
  if (i%%2==0) {family1[i]=""}
  if (i%%2==1) {family2[i]=""}
}
N=ifelse(estnetwork1==0,0,1)
g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "directed",diag = F)
V(g2)$degree = degree(g2)
V(g2)$name=family
node_list=get.data.frame(g2,what="vertices")
edge_list=get.data.frame(g2,what="edges")

plot_data=data.frame()
for (i in rownames(estnetwork1)) {
  for (j in colnames(estnetwork1)) {
    a=data.frame(from=i,to=j,coefficient=abs(estnetwork1[i,j]))
    plot_data=rbind(plot_data,a)
  }
}








#png(filename = paste("networkcommon_he", ".png", sep = ""))[1,]
ggplot(plot_data, aes(x = from, y = to, fill = coefficient))  +
  geom_tile() +
  theme_bw(base_size = 8, base_family = "serif", base_line_size = 1) +
  scale_x_discrete(labels=family2) +
  scale_y_discrete(labels=family1) +
  labs(x="Taxon",y="Taxon")+
  theme(
    # Rotate the x-axis labels so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0),
    # Force the plot into a square aspect ratio
    aspect.ratio = 1,
    # Hide the legend (optional)
    #legend.position = "none"
    )
  theme(legend.title = element_text(colour="white", size = 16, face='bold'))
#dev.off()

  ##
  cc=cluster_louvain(g2, weights = NULL, resolution = 1)
  ll=c()
  for (i in 1:length(cc)) {
    ll=c(ll,length(cc[[i]]))
  }
  ll=sapply(cc, length)
  index=which(ll>1)
  commu=list()
  for (i in 1:length(index)) {
    commu[[i]]=family[cc[[index[i]]]]
  }


  belong=c()
  for (i in family) {
    a=ifelse( i%in% commu[[1]],1,ifelse(i %in% commu[[2]], 2,
                                        ifelse(i %in% commu[[3]],3,
                                               ifelse(i %in% commu[[4]], 4, ifelse( i%in% commu[[5]], 5,0)))) )
    belong=c(belong,a)
  }

