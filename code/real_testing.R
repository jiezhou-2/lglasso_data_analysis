data=read.csv(file = "final.csv", header = T)[,-1]
subject=data[,1]
age=data[,2]
agelevel=levels(factor(age))
rabundance=data[,3:ncol(data)]
lrdata=matrix(nrow = nrow(rabundance), ncol = ncol(rabundance)-1)
for (j in 1:ncol(lrdata)) {
  index=which(rabundance[,j]!=0)
  impute=min(rabundance[index,j])/10
  a1=ifelse(rabundance[,j]!=0,rabundance[,j],impute)
  a2=rabundance[,ncol(rabundance)]
  lrdata[,j]=log(a1/a2)
}
datalist=vector("list", 17)
for (k in 1:(length(datalist))){
  index=which(age==agelevel[k])
  D=lrdata[index,]
  ave=apply(lrdata[index,],2,mean)
  datalist[[k]]=t(apply(D,1,"-",ave))
}
a=rep(0,16)
for (j in 1:16) {
a[j]=twosample(datalist[[j]],datalist[[j+1]])
}

pv=1-pnorm(abs(a))
divide=c(1,1,2,2,2,2,3,3,4,4,4,4,4,4,4,4,4)
