getResidMatrix=function(ymat, X=NULL, mask){
  mask.index=which(mask!=0)
  if (X=NULL){
    for (j in mask.index){
      ymat[,j]=ymat[,j]/sum(ymat[,j]^2,na.rm=T)
    }
  }
  return(ymat)
}