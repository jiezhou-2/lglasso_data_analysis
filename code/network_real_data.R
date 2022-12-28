library("devtools")
load_all("../../lglasso")
data=read.csv(file = "C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso_data_analysis/data/real/ddata.csv")[,-1]
bic1=c()
rho1=seq(1,13,length=5)
Te=1
s=cov(data[,-c(1,2)])
network1=vector("list",length = length(rho1))
for (j in 1:length(rho1)){
  print(j)
result1=lglasso(data = data, rho = rho1[j],heter=F)
network1[[j]]=result1$omega
bb1=-2*result1$ll +0.5*length(which(result1$omega!=0))*log(nrow(data))+0.5* length(which(result1$omega!=0))*log(ncol(data)-2)/Te
bic1=c(bic1,bb1)
}
save(network1,file="network1.rd")
save(bic1,file="bic1.rd")



















