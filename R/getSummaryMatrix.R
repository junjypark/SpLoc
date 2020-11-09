getResidMatrix=function(ymat, X=NULL, mask){
  mask.index=which(mask!=0)
  if (is.null(X)){
    for (j in mask.index){
      ymat[,j]=ymat[,j]/sum(ymat[,j]^2,na.rm=T)
    }
  } else{
    X=as.matrix(X)
    n=nrow(X); p=ncol(X)
    P=tcrossprod(tcrossprod(X,solve(crossprod(X))),X)
    Q=1-apply(P,1,sum)
    ymat=Q*ymat
  }
  ymat[is.nan(ymat)]=0
  return(ymat)
}

