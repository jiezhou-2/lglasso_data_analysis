bicfunction=function(data,G,nu=1){
  if (det(G)<=0){
    warning("G is not positive definite matrix!")
    return(bic=NA)
  }else{
    n=nrow(data)
    k=length(which(G!=0))/2
    yy=scale(as.matrix(data[,-c(1,2)]))
    ll=-0.5*sum(apply(yy,1,function(x)  t(x) %*% G %*% x ))+0.5*n*log(det(G))
    bic=-2*ll+k*log(n)*nu
  }
  return(bic=bic)
}
