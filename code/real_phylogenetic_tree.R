##smalltree contains the evolution information
taxa=read.csv(file="taxa.csv",header = T)[,-1]
taxa[which(taxa=="Escherichia/Shigella",arr.ind = T)]="Escherichia.Shigella"
#target=c(commu1,commu2,commu3,commu4,commu5)
family[2]="Escherichia.Shigella"
target=family#here data is generated from the file realdata.R
smalltree=data.frame()
for (i in 1:length(target)) {
  index=which(taxa[,6]==target[i])
  print(length(index))
  a=taxa[index,]
  smalltree=rbind(smalltree,a)
}
smalltree=unique(smalltree)



##generate file for drawing the phylogenetic  tree used in software graphlan
graphlanTree=data.frame(matrix(nrow = nrow(smalltree),ncol = 1))
for (i in 1:nrow(smalltree)){
  index=which(is.na(smalltree[i,]))
  if(length(index)==0){
    pure=smalltree[i,]
  }else{
    pure=smalltree[i,-index]
  }
  a=pure[1]
  for (j in 2:length(pure)) {
    a=paste(a,pure[j],sep = ".")
  }
  graphlanTree[i,]=a
}
graphlanTree=unique(graphlanTree)
write.csv(graphlanTree,file="graphlanTree.csv")


##generate edge from a named path 
gedge=function(pure){
  n=length(pure)
  edgematrix=as.data.frame(matrix(nrow = n-1, ncol = 2))
  for (i in 1:(n-1)){
    edgematrix[i,]=c(pure[i],pure[i+1])
  }
  return(edgematrix)
}
##generate a graph using evolution information
edgeset=data.frame()
for (k in 1:nrow(smalltree)){
  index=which(is.na(smalltree[k,]))
  if(length(index)==0){
    pure=smalltree[k,]
  }else{
    pure=smalltree[k,-index]
  }
  edgeset=rbind(edgeset,gedge(pure))
}
treegraph=graph_from_edgelist(as.matrix(edgeset),directed = F)
M=as_adjacency_matrix(treegraph)
