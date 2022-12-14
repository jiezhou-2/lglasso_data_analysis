
sim_heter=function(p,prob,alpha,age){
  print("heter data are generated")
  m=length(age)
  data=vector("list",m)
  tau=rexp(length(alpha),rate = alpha)
  K=BDgraph::bdgraph.sim(p=p,n=1,type="Gaussian",prob=prob)$K
  Sigma=solve(K)
  sqK=chol(Sigma)
  for (i in 1:m) {
    n=length(age[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK)%*%error1)
    a[1,]=error2[1,]
  for (t in 2:n) {
    coe=exp(-tau[i]*abs(age[[i]][t]-age[[i]][t-1]))
    a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
  }
    data[[i]]=cbind(i,age[[i]],a)
  }
  result=list(data=data,tau=tau,alpha=alpha,precision=K)
  return(result)
}





sim_homo=function(p,prob,tau,age){
  print("homo data are generated")
  m=length(age)
  data=vector("list",m)
  K=BDgraph::bdgraph.sim(p=p,n=1,type="Gaussian",prob=prob)$K
  Sigma=solve(K)
  sqK=chol(Sigma)
  for (i in 1:m) {
    n=length(age[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK)%*%error1)
    a[1,]=error2[1,]
    for (t in 2:n) {
      coe=exp(-tau*abs(age[[i]][t]-age[[i]][t-1]))
      a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
    }
    data[[i]]=cbind(i,age[[i]],a)
  }
  result=list(data=data,tau=tau,precision=K)
  return(result)
}



comparison=function(real, estimate){
  real=real+t(real)
  diag(real)=1
  estimate=estimate+t(estimate)
  diag(estimate)=1
  N1=ifelse(abs(real)<=10^(-5),0,1)
  N2=ifelse(abs(estimate)<=10^(-5),0,1)
  if (any(dim(N1)!=dim(N2)))
    stop("Two matrixes should have the same dimension")
  p=nrow(real)
  real=(sum(N1)-p)/2
  null=p*(p-1)/2-real
  select=(sum(N2)-p)/2
  real_select=(sum(N2[N1==1])-p)/2
  fause_select=sum(N2[N1==0])/2
  if (real==0){
    stop("Real network has no edge. Please regenerate the real network")
    }
  if (select==0) {
    stop("Estimated network has no edge. Please decrease the tuning parameter.")
    }
  TPR=real_select/real
  FPR=fause_select/null
  #FPR=fause_select/select
  aa=as.numeric(c(TPR, FPR))
  return(aa)
}





addition=function(data,lambda){
  data=data[,-c(1,2)]
  p=ncol(data)
  n=nrow(data)
  arr=matrix(0,nrow=p,ncol=p)

  for (i in 1:p) {
    x=as.matrix(data[,-i])
    y=data[,i]
    result=glmnet::glmnet(x=x,y=y, family="gaussian",lambda = lambda)
    bb=as.matrix(result$beta)
    arr[i,-i]=bb
  }
  web=arr+t(arr)
  return(web)
}




power_compare1=function(m,n,p,coe,l,rho,prob,heter,community2=F,uu=c(0,0),zirate=c(0.2,0)){
  results=vector("list",length = 5) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho[[1]]))
  results[[2]]=vector("list",length = length(rho[[2]]))
  results[[3]]=vector("list",length = length(rho[[3]]))
  results[[4]]=vector("list",length = length(rho[[4]]))
  results[[5]]=vector("list",length = length(rho[[5]]))
  RR=vector("list",length = 5)
  RR[[1]]=matrix(nrow=l,ncol = 2)
  RR[[2]]=matrix(nrow=l,ncol = 2)
  RR[[3]]=matrix(nrow=l,ncol = 2)
  RR[[4]]=matrix(nrow=l,ncol = 2)
  RR[[5]]=matrix(nrow=l,ncol = 2)

  age=vector("list",m) #generate the time points container
  for (i in 1:5) {
    for (j in 1:length(rho[[i]])) {
      results[[i]][[j]]= matrix(nrow = l,ncol = 2)
    }
  }


  for (h in 1:l){
    print(paste0("The ",h,"th"," simulation:"))
    # subject level covariates
    x1=sample(x=c(0,1),size=m,prob=c(0.5,0.5),replace = T)
    x2=runif(m,min = 0,max = 1)
    alpha=exp(coe[1]+coe[2]*x1+coe[3]*x2)
    x=cbind(x1,x2)
    ## generate the observation time for each subjects
    for (k in 1:m) {
      a1=rpois(n=n,lambda = 1) # the space between observations
      age[[k]][1]=max(a1[1],0.5)
      for (i in 2:n) {
        age[[k]][i]=age[[k]][i-1]+max(a1[i-1],0.5)
      }
    }

    ## generate the network data

    if (heter==TRUE & community2==F){
      ss=sim_heter(p = p,prob=prob,alpha = alpha,age = age,zirate=zirate)
    }
    if (heter==TRUE & community2==T){
      ss=sim_2heter(p=p,prob=prob,alpha1=exp(uu[1]),alpha2=exp(uu[2]),age=age,zirate = zirate)
    }
    if (heter==F & community2==F){
      ss=sim_homo(p = p,prob=prob,tau = 1/alpha[1],age = age,zirate=zirate)
    }
    if (heter==F & community2==T){
      ss=sim_2homo(p=p,prob=prob,tau1=1/exp(uu[1]),tau2=1/exp(uu[2]),age=age,zirate = zirate)
    }

    simdata=ss$data
    lower=0.01
    upper=20
    graph=ss$precision
    dd=do.call(rbind,simdata)
    id=unique(dd[,1])
    covariate=cbind(id,x)
    Tem=log(p)/(2*log(1/prob-1))

    bic1=c()
    bic2=c()
    bic3=c()
    bic4=c()
    bic5=c()
    aa=vector("list",length(rho[[1]]))
    aa1=vector("list",length(rho[[2]]))
    aa2=vector("list",length(rho[[3]]))
    aa3=vector("list",length(rho[[4]]))
    aa4=vector("list",length(rho[[5]]))

    ## lglasso
    for (j in 1:length(rho[[1]])) {
      if (heter==TRUE){
        if (length(which(coe==0))==2){
          aa[[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j],heter=T)
        }else{
          aa[[j]]=lglasso(data = dd,x=covariate, rho = 0.5*rho[[1]][j],heter=T)
        }
      }else{
        aa[[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j])
      }
      bb1=-2*aa[[j]]$ll +  0.5*length(which(aa[[j]]$omega!=0))*log(nrow(dd))
      +0.5*length(which(aa[[j]]$omega!=0))*log(p)/Tem
      bic1=c(bic1,bb1)
      results[[1]][[j]][h,]=as.numeric(comparison(graph,aa[[j]]$omega))
    }

    print(paste0("lglasso's bic is ",bic1))

    ## glasso

    for (j in 1:length(rho[[2]])) {
      ##estiamte the network based on glasso
      s=cov(dd[,-c(1,2)])
      aa1[[j]]=glasso(s=s,rho=rho[[2]][j])$wi
      bb2=bicfunction(data=dd,G=as.matrix(aa1[[j]]),Tem=Tem)
      bic2=c(bic2,bb2)
      results[[2]][[j]][h,]=as.numeric(comparison(graph,aa1[[j]]))
    }


    ## nh

    for (j in 1:length(rho[[3]])) {
      bb2=addition(data=dd,lambda=rho[[3]][j])
      aa2[[j]]=mle_net(data=dd,priori=bb2)
      bb3=bicfunction(data=dd,G=aa2[[j]],Tem=Tem)
      bic3=c(bic3,bb3)
      results[[3]][[j]][h,]=as.numeric(comparison(graph,aa2[[j]]))
    }
    ## CO1
    for (j in 1:length(rho[[4]])) {
      aa3[[j]]=selectFast(s,family="C01",K=2*rho[[4]][j])$C01$G
      bb4=bicfunction(data=dd,G=aa3[[j]],Tem=Tem)
      bic4=c(bic4,bb4)
      results[[4]][[j]][h,]=as.numeric(comparison(graph,aa3[[j]]))

    }

    ## la
    for (j in 1:length(rho[[5]])) {
      aa4[[j]]=selectFast(s,family="LA",K=2*rho[[5]][j])$LA$G
      bb5=bicfunction(data=dd,G=aa4[[j]],Tem=Tem)
      bic5=c(bic5,bb5)
      results[[5]][[j]][h,]=as.numeric(comparison(graph,aa4[[j]]))
    }

    G1=aa[[which.min(bic1)]]$omega
    #print(bic2)
    G2=aa1[[which.min(bic2)]]
    #print(bic3)
    G3=aa2[[which.min(bic3)]]
    G4=aa3[[5]]
    G5=aa4[[5]]


    RR[[1]][h,]=as.numeric(comparison(graph,G1))
    RR[[2]][h,]=as.numeric(comparison(graph,G2))
    RR[[3]][h,]=as.numeric(comparison(graph,G3))
    RR[[4]][h,]=as.numeric(comparison(graph,G4))
    RR[[5]][h,]=as.numeric(comparison(graph,G5))
  }


  bb=vector("list",length = 5)
  bb[[1]]=matrix(nrow = length(rho[[1]]),ncol=2)
  bb[[2]]=matrix(nrow = length(rho[[2]]),ncol=2)
  bb[[3]]=matrix(nrow = length(rho[[3]]),ncol=2)
  bb[[4]]=matrix(nrow = length(rho[[4]]),ncol=2)
  bb[[5]]=matrix(nrow = length(rho[[5]]),ncol=2)
  names(bb)=c("lglasso","glasso","nh","co","la")
  for (i in 1:5) {
    for (j in 1:length(rho[[i]])) {
      iimatrix=results[[i]][[j]]
      index1=which( is.na(iimatrix[,1]))
      index2=which(is.na(iimatrix[,2]))
      index=union(index1,index2)
      if (length(index)==l){
        warning(paste0(" All results are NA for rho= ", i,"+",j,rho[[i]][j]))
        next
      }
      ff=if(length(index)==0){
        ff=iimatrix
      }else{
        ff=iimatrix[-index,]
      }
      bb[[i]][j,]=apply(ff, 2, mean)
    }
  }

  rr=matrix(nrow=5,ncol=2)

  for (i in 1:5) {
    rr[i,]=apply(RR[[i]], 2, mean)
  }


  return(list(whol=bb,selec=rr,resul=results))
}




