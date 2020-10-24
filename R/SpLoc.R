SpLoc=function(NNmatrix, ymat, nperm=5000, alpha=0.05, seed=NULL, is.sparse=F){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (is.null(seed)){
    stop("Specifying a seed value is required.")
  }
  if (is.sparse){
    index=NULL
    for (j in 1:ncol(NNmatrix)){
      if(sum(NNmatrix[,j])!=0){ index=c(index,j)}
    }
    NNmatrix=NNmatrix[,index]
    ymat=ymat[index,]
  }

  pU=big.matrix(nrow(NNmatrix), nperm, type = "double")
  out=SpLocC(NNmatrix, ymat, nperm, alpha, seed, pU@address)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
  out$seed=seed
  out$nperm=nperm
  return(out)
}

combine1=function(lst, alpha=0.05){
  n=length(lst)
  seed=do.call("c",lapply(lst, function(x){x$seed}))
  nperm=do.call("c",lapply(lst, function(x){x$nperm}))
  
  if (length(unique(seed))>1){stop("Use the same seed for every Sploc output.")}
  if (length(unique(nperm))>1){stop("Use the same number of permutations for every SpLoc output.")}
  
  Tstat=do.call("c",lapply(lst, function(x){x$Tstat}))
  permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)
  threshold=quantile(permMax, 1-alpha)
  pvalue=(1+sum(c(permMax)>max(Tstat,na.rm=T)))/(1+nperm[1])
  
  return(list(
    threshold=threshold,
    Tstat=Tstat,
    permMax=permMax,
    pvalue=pvalue,
    seed=seed[1],
    nperm=nperm[1]
  ))
}
