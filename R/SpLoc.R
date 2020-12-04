SpLoc=function(ymat, NNmatrix=NULL, group=NULL, nperm=1000, alpha=0.05, seed=NULL, 
               is.sparse=F,partition=F, npartition=1, parallel=F, ncores=1){
  if (is.null(NNmatrix)){
    print("As NNmatrix is not specified, SpLoc conducts massive univariate analysis.")
    if (is.null(group)){
      return(MassiveMean(ymat=ymat,
                         nperm=nperm, 
                         alpha=alpha, 
                         seed=seed,
                         partition=partition, 
                         npartition = npartition,
                         parallel=parallel,
                         ncores=ncores)
                         )
    } else{
      return(MassiveDiff(ymat=ymat,
                         group=group,
                         nperm=nperm, 
                         alpha=alpha, 
                         seed=seed,
                         partition=partition, 
                         npartition = npartition,
                         parallel=parallel,
                         ncores=ncores)
                         ) 
    }
    
  } else{
    if (is.null(group)){
      return(
        SpLocMean(ymat=ymat, 
                  NNmatrix=NNmatrix, 
                  nperm=nperm, 
                  alpha=alpha, 
                  seed=seed,
                  is.sparse = is.sparse, 
                  partition=partition, 
                  npartition = npartition,
                  parallel=parallel,
                  ncores=ncores)
      )
    }
    else{
      return(
        SpLocDiff(ymat=ymat, 
                  NNmatrix=NNmatrix, 
                  group=group, 
                  nperm=nperm, 
                  alpha=alpha, 
                  seed=seed,
                  is.sparse = is.sparse, 
                  partition=partition, 
                  npartition = npartition,
                  parallel=parallel, 
                  ncores=ncores)
      )
    }
  }
}


SpLocMean=function(ymat, NNmatrix, nperm=1000, alpha=0.05, seed=NULL, 
                   is.sparse=F,partition=F, npartition=1, parallel=F, ncores=1){
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
        SpLocMeanC(ymat, NNList[[i]], nperm, alpha, seed)
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=SpLocMeanC(ymat, NNList[[i]], nperm, alpha, seed)
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
    
  } else{
    out=SpLocMeanC(ymat, NNmatrix, nperm, alpha, seed)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    return(out)    
  }
}


SpLocDiff=function(ymat, NNmatrix, group, nperm=1000, alpha=0.05, seed=NULL, 
                   is.sparse=F,partition=F, npartition=1, parallel=F, ncores=1){
  if (!all.equal(sort(unique(group)),c(-1,1))){
    stop("group should have either 1 or -1.")
  }
  if (length(group)!=ncol(ymat)){
    stop("The number of elements in group does not match with the number of columns in ymat.")
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
        SpLocDiffC(ymat, NNList[[i]], group, nperm, alpha, seed)
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=SpLocDiffC(ymat, NNList[[i]], group, nperm, alpha, seed)
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
    
  } else{
    out=SpLocDiffC(ymat, NNmatrix, group, nperm, alpha, seed)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    return(out)    
  }
}


MassiveMean=function(ymat, nperm=1000, alpha=0.05, seed=NULL, 
                   partition=F, npartition=1, parallel=F, ncores=1){
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (is.null(seed)){
    stop("Specifying a seed value is required.")
  }

  if (isTRUE(partition)){
    ymatList=list()
    len=ceiling(nrow(ymat)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(ymat))
      ymatList[[i]]=ymat[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("SpLoc"),.noexport = "SpLocC" )%dopar%{
        MassiveMeanC(ymatList[[i]], nperm, alpha, seed)
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=MassiveMeanC(ymatList[[i]], nperm, alpha, seed)
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
    
  } else{
    out=MassiveMeanC(ymat, nperm, alpha, seed)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
  }
  
  return(out)    
}


MassiveDiff=function(ymat, group, nperm=1000, alpha=0.05, seed=NULL, 
                   is.sparse=F,partition=F, npartition=1, parallel=F, ncores=1){
  if (!all.equal(sort(unique(group)),c(-1,1))){
    stop("group should have either 1 or -1.")
  }
  if (length(group)!=ncol(ymat)){
    stop("The number of elements in group does not match with the number of columns in ymat.")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (is.null(seed)){
    stop("Specifying a seed value is required.")
  }
  
  out=MassiveDiffC(ymat, group, nperm, alpha, seed)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
  out$seed=seed
  return(out)    
}

