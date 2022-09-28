data=read.csv(file = "ddata.csv")[,-1]
alpha0=0.1
tau0=10
p=ncol(data)-2
omega0=diag(p)
rho=60
#rho_homo=70
#For a
#rho=exp(0.1)
tole=0.1
#a=temglasso_heter(data = data, alpha0 = alpha0,omega0 = omega0,rho = rho,error = error,tole=tole,type="abs")
a=longraph(data = data,alpha0 = alpha0,omega0=omega0,type="abs", rho = 48,tole=0.01, lower = 0.1, upper = 7)
#path analysis
rho1=seq(0.1,3.252632,length=10)
#rho=seq(45,55,length=40)
rho2=seq(3.3,45,length=20)
rho3=seq(45,55,length=40)
rho4=seq(55,70,length=20)
rho=c(rho1,rho2,rho3,rho4)
pathlist=vector("list",length(rho))
for (i in 1:length(rho)) {
  print(i)
  a=longraph(data = data,alpha0 = alpha0,omega0=omega0,type="abs", rho = rho[i],tole=0.01, lower = 0.1, upper = 7)
  N=ifelse(a$omega==0,0,1)
  pathlist[[i]]=N
}
u=c()
for (i in 1:length(pathlist)) {
  u=c(u,(sum(pathlist[[i]])-82)/2)
}



##for heterogeneous model
mlenet_he=mle(data =d, priori = a$omega)
alpha_tau_mle=mle_alpha(data = data,alpha0 = 0.5,omega=mlenet_he, type="abs", tole=0.01, lower = 0.2, upper = 7)
tau_he=alpha_tau_mle$tau
dottau=data.frame(x="Cohort",y=tau_he)
theme_set(theme_classic()+
            theme(legend.position = "top"))
e=ggplot(dottau,aes(x=y))
e+geom_histogram(aes(y=stat(density)),binwidth=0.3, colour="black",fill="white")+
  geom_density(alpha=0.2,fill="#E7B800")+
  labs(x="Autocorrelation Parameter",y="Density")

e=ggplot(dottau,aes(x=x,y=y))
e+geom_boxplot(width=0.3)+
  geom_dotplot(binaxis = 'y',stackdir = 'center')+
  xlim("Cohort")+
labs(x="",y="Autocorrelation parameter")

hist(tau_he,)
title("")
plot(tau_he)







##
family=colnames(data)[-c(1,2)]
family1=family
l=length(family)
for (i in 1:l) {
  if (i%%2==0) {family1[i]=""} 
}
N=ifelse(a$omega==0,0,1)
N=pathlist[[23]]
write.csv(family,file="genus.csv")
write.csv(N,file="adjacency_matrix.csv")
g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "directed",diag = F)
V(g2)$degree = degree(g2)
V(g2)$name=family
node_list=get.data.frame(g2,what="vertices")
edge_list=get.data.frame(g2,what="edges")
group=rep(1,nrow(edge_list))
edge_list=cbind(edge_list,group)
all_nodes <- sort(node_list$name)
plot_data <- edge_list %>% mutate(
  to = factor(to, levels = all_nodes),
  from = factor(from, levels = all_nodes))
png(filename = paste("networkcommon_he", ".png", sep = ""))[1,]
ggplot(plot_data, aes(x = from, y = to, fill = group))  +
  geom_raster() +
  theme_bw(base_size = 11, base_family = "serif", base_line_size = 1) +
  # Because we need the x and y axis to display every node,
  # not just the nodes that have connections to each other,
  # make sure that ggplot does not drop unused factor levels
  # scale_x_discrete(drop = FALSE) +
  scale_x_discrete(labels=family1) +
  scale_y_discrete(labels=family1) +
  labs(x="Taxon",y="Taxon")+
  theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0),
    # Force the plot into a square aspect ratio
    aspect.ratio = 1,
    # Hide the legend (optional)
    legend.position = "none")
dev.off()

###transform the network to community
N=as.matrix(read.csv(file = "adjacency_matrix.csv")[,-1])
rownames(N)=family
colnames(family)
graphHeter=graph_from_adjacency_matrix(N,mode = "undirected")
#amatrix=as_adjacency_matrix(graphHeter) 
community=cluster_louvain(graph = graphHeter)
member=communities(community)
commu1=family[which(community$membership==6)]
commu2=family[which(community$membership==7)]
commu3=family[which(community$membership==8)]
commu4=family[which(community$membership==12)]
commu5=family[which(community$membership==13)]
commu6=family[which(community$membership==16)]
commu7=family[which(community$membership==21)]
commu8=family[which(community$membership==22)]
commu=c(commu1,commu2,commu3,commu4,commu5,commu6,commu7,commu8)
annot_1=data.frame(matrix(nrow = 2*length(commu),ncol = 3))
for (i in 1:length(commu)) {
  annote1[2*i-1,1]=commu[i]
  annote1[2*i-1,2]="clade_marker_color"
  if (commu[i]%in%commu1){annote1[2*i-1,3]="b"}
  if (commu[i]%in%commu2){annote1[2*i-1,3]="r"}
  if (commu[i]%in%commu3){annote1[2*i-1,3]="#FF4433"}
  if (commu[i]%in%commu4){annote1[2*i-1,3]=" #00AAAA"}
  if (commu[i]%in%commu5){annote1[2*i-1,3]="k"}
  if (commu[i]%in%commu6){annote1[2*i-1,3]="#008080"}
  if (commu[i]%in%commu7){annote1[2*i-1,3]="#f1cbff"}
  if (commu[i]%in%commu8){annote1[2*i-1,3]="#c99789"}
  annote1[2*i,1]=commu[i]
  annote1[2*i,2]="clade_marker_size"
  annote1[2*i,3]=20
}
write.csv(annote1,file = "annot_1.csv")


a_g=glasso::glasso(s=cov(data[,-c(1,2)]),rho = 10)$wi 
family=colnames(data)[-c(1,2)]
N=ifelse(a_g==0,0,1)
g2=graph_from_adjacency_matrix(adjmatrix = t(N), mode = "directed",diag = F)
V(g2)$degree = degree(g2)
V(g2)$name=family
node_list=get.data.frame(g2,what="vertices")
edge_list=get.data.frame(g2,what="edges")
group=rep(1,nrow(edge_list))
edge_list=cbind(edge_list,group)
all_nodes <- sort(node_list$name)
plot_data <- edge_list %>% mutate(
  to = factor(to, levels = all_nodes),
  from = factor(from, levels = all_nodes))
png(filename = paste("networkcommon_g", ".png", sep = ""))
ggplot(plot_data, aes(x = from, y = to, fill = group,xlab="Microbe",ylab="Microbes"))  +
  geom_raster() +
  theme_bw() +
  # Because we need the x and y axis to display every node,
  # not just the nodes that have connections to each other,
  # make sure that ggplot does not drop unused factor levels
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x="Microbe",y="Microbe")+
  theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0),
    # Force the plot into a square aspect ratio
    aspect.ratio = 1,
    # Hide the legend (optional)
    legend.position = "none")
dev.off()

