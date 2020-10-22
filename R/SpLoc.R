SpLoc=function(NNmatrix, ymat, nperm=5000, alpha=0.05, seed=NULL){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it")
  }

  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels)")
  }

  if (is.null(seed)){
    print("Specifying a seed value is recommended combining multiple NN matrix.")
  } else{
    set.seed(seed)
  }

  out=SpLocC(NNmatrix, ymat, nperm, alpha)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat)))/(1+nperm)

  return(out)
}

ClusterSearch=function(tstat, thres, NNmatrix){
  if (length(tstat)!=nrow(NNmatrix)){
    stop("The number of rows in NNmat needs to be the same as the length of tstat,")
  }

  sig.ind=which(tstat>thres)
  tstat=tstat[sig.ind]
  NNmatrix=NNmatrix[sig.ind,,drop=F]

  if (nrow(NNmatrix)==0){bool=F}else{bool=T}

  clust=1
  sig=NULL
  while(bool){
    ind=which(tstat==max(tstat))[1]
    sig.vertices=which(NNmatrix[ind,]!=0)
    sig=c(sig, sig.vertices)
    out.set=which(apply(NNmatrix[,sig.vertices,drop=F],1,max)>0)

    tstat=tstat[-out.set]
    NNmatrix=NNmatrix[-out.set,,drop=F]

    print(paste0("Cluster ", clust, " with ", length(sig.vertices), " voxel(s) is selected (total=", length(sig), ")" ))
    clust=clust+1

    if (length(tstat)==0){bool=F}
  }

  return(sig)
}

