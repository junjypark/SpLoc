SpLoc=function(NN, ymat, nperm=5000, alpha=0.05, seed=NULL){
  if (length(which(is(NN)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it")
  }

  if ( ncol(NN)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels)")
  }

  if (is.null(seed)){
    print("Specifying a seed value is recommended combining multiple NN matrix.")
  } else{
    set.seed(seed)
  }
  print(dim(NN))
  out=SpLocC(NN, ymat, nperm, alpha)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat)))/(1+nperm)

  return(out)
}

ClusterSearch=function(tstat, thres, NNmat){
  if (length(tstat)!=nrow(NNmat)){
    stop("The number of rows in NNmat needs to be the same as the length of tstat,")
  }

  sig.ind=which(tstat>thres)
  tstat=tstat[sig.ind]
  NNmat=NNmat[sig.ind,,drop=F]

  if (nrow(NNmat)==0){bool=F}else{bool=T}

  clust=1
  sig=NULL
  while(bool){
    ind=which(tstat==max(tstat))[1]
    sig.vertices=which(NNmat[ind,]!=0)
    sig=c(sig, sig.vertices)
    out.set=which(apply(NNmat[,sig.vertices,drop=F],1,max)>0)

    tstat=tstat[-out.set]
    NNmat=NNmat[-out.set,,drop=F]

    print(paste0("Cluster ", clust, " with ", length(sig.vertices), " voxel(s) is selected (total=", length(sig), ")" ))
    clust=clust+1

    if (length(tstat)==0){bool=F}
  }

  return(sig)
}

