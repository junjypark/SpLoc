SpLoc=function(NNmatrix, ymat, group=NULL, nperm=1000, alpha=0.05, seed=NULL, 
               is.sparse=F,partition=F, npartition=NULL, parallel=F, ncores=1){
  if (is.null(group)){
    return(SpLocMean(...))
  }
  else{
    return(SpLocDiff(...))
  }
}

SpLocMean=function(NNmatrix, ymat, nperm=1000, alpha=0.05, seed=NULL, 
                is.sparse=F,partition=F, npartition=NULL, parallel=F, ncores=1){
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
  if (isTRUE(partition)){
    if (is.null(npartition)){
      partition=FALSE
      print("partition set to be FALSE as npartition is not specified.")
    }
    else if (npartition==1){
      partition=FALSE
      print("partition set to be FALSE as npartition is 1.")
    }
  }
  
  if (isTRUE(partition)){
    NNList=list()
    len=ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(NNmatrix))
      NNList[[i]]=NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("SpLoc"),.noexport = "SpLocC" )%dopar%{
        pU=big.matrix(nrow(NNList[[i]]), nperm, type = "double")
        pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
        SpLocMeanC(NNList[[i]], ymat, nperm, alpha, seed, pU@address, pY@address)
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        pU=big.matrix(nrow(NNList[[i]]), nperm, type = "double")
        pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
        result[[i]]=SpLocMeanC(NNList[[i]], ymat, nperm, alpha, seed, pU@address,pY@address)
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
    
  } else{
    pU=big.matrix(nrow(NNmatrix), nperm, type = "double")
    pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
    out=SpLocMeanC(NNmatrix, ymat, nperm, alpha, seed, pU@address, pY@address)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    return(out)    
  }
}


SpLocDiff=function(NNmatrix, ymat, group, nperm=1000, alpha=0.05, seed=NULL, 
               is.sparse=F,partition=F, npartition=NULL, parallel=F, ncores=1){
  if (!all.equal(sort(unique(group)),c(-1,1))){
    stop("group should have either 1 or -1")
  }
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
  if (isTRUE(partition)){
    if (is.null(npartition)){
      partition=FALSE
      print("partition set to be FALSE as npartition is not specified.")
    }
    else if (npartition==1){
      partition=FALSE
      print("partition set to be FALSE as npartition is 1.")
    }
  }
  
  
  if (isTRUE(partition)){
    NNList=list()
    len=ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(NNmatrix))
      NNList[[i]]=NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("SpLoc"),.noexport = "SpLocC" )%dopar%{
        pU=big.matrix(nrow(NNList[[i]]), nperm, type = "double")
        pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
        SpLocDiffC(NNList[[i]], ymat, group, nperm, alpha, seed, pU@address, pY@address)
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        pU=big.matrix(nrow(NNList[[i]]), nperm, type = "double")
        pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
        result[[i]]=SpLocDiffC(NNList[[i]], ymat, group, nperm, alpha, seed, pU@address, pY@address)
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
    
  } else{
    pU=big.matrix(nrow(NNmatrix), nperm, type = "double")
    pY=big.matrix(ncol(NNmatrix), nperm, type = "double")
    out=SpLocDiffC(NNmatrix, ymat, group, nperm, alpha, seed, pU@address, pY@address)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    return(out)    
  }
}
