getResidMatrix=function(ymat, X=NULL, mask){
  mask.index=which(mask!=0)
  if (is.null(X)){
    for (j in mask.index){
      ymat[,j]=ymat[,j]/sum(ymat[,j]^2,na.rm=T)
    }
  }
  ymat[is.nan(ymat)]=0
  return(ymat)
}