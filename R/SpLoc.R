SpLoc=function(ymat, NNmatrix=NULL, group=NULL, nperm=10000, alpha=0.05, alternative=c("two.sided","less", "greater"), seed=NULL, 
               partition=T, npartition=NULL, parallel=F, ncores=1){
  if (is.null(NNmatrix)){
    print("As NNmatrix is not specified, SpLoc conducts massive univariate analysis.")
    if (is.null(group)){
      fit=MassiveMean(ymat=ymat,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit)
    } else{
      fit=MassiveDiff(ymat=ymat,
                      group=group,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit) 
    }
  } else{
    if (is.null(group)){
      fit=SpLocMean(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel,
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
    else{
      fit=SpLocDiff(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    group=group, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel, 
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
  }
}


SpLocMean=function(ymat, NNmatrix, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                   partition=F, npartition=1, parallel=F, ncores=1){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    print("Conducting the two-sided test as alternative has not been specified...")
    }
  if (alternative=="two.sided"){ side=2 }
  if (alternative=="greater"){ side=1 }
  if (alternative=="less"){ side=-1 }
  if (is.null(seed)){ seed=sample(1e6,1) }
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    else if (npartition==1){
      partition=FALSE
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
        fit=SpLocMeanC(ymat, NNList[[i]], nperm, seed, side)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=SpLocMeanC(ymat, NNList[[i]], nperm, seed, side)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
  } else{
    out=SpLocMeanC(ymat, NNmatrix, nperm, seed, side)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    out$alternative=alternative
    return(out)    
  }
}


SpLocDiff=function(ymat, NNmatrix, group, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                   partition=F, npartition=1, parallel=F, ncores=1){
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
  if (length(alternative)>1){ 
    alternative="two.sided"
    print("Conducting the two-sided test as alternative has not been specified...")
  }
  if (alternative=="two.sided"){ side=2 }
  if (alternative=="greater"){ side=1 }
  if (alternative=="less"){ side=-1 }
  if (is.null(seed)){ seed=sample(1e6,1) }
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    else if (npartition==1){
      partition=FALSE
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
        fit=SpLocDiffC(ymat, NNList[[i]], group, nperm, seed, side)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=SpLocDiffC(ymat, NNList[[i]], group, nperm, seed, side)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
  } else{
    out=SpLocDiffC(ymat, NNmatrix, group, nperm, seed, side)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    out$alternative=alternative
    return(out)    
  }
}


MassiveMean=function(ymat, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                   partition=F, npartition=1, parallel=F, ncores=1){
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    print("Conducting the two-sided test as alternative has not been specified...")
  }
  if (alternative=="two.sided"){ side=2 }
  if (alternative=="greater"){ side=1 }
  if (alternative=="less"){ side=-1 }
  if (is.null(seed)){ seed=sample(1e6,1) }
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    
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
        fit=MassiveMeanC(ymatList[[i]], nperm, seed, side)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=MassiveMeanC(ymatList[[i]], nperm, seed, side)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
  } else{
    out=MassiveMeanC(ymat, nperm, seed, side)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    out$alternative=alternative
  }
  
  return(out)    
}


MassiveDiff=function(ymat, group, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                     partition=F, npartition=1, parallel=F, ncores=1){
  if (!all.equal(sort(unique(group)),c(-1,1))){
    stop("group should have either 1 or -1.")
  }
  if (length(group)!=ncol(ymat)){
    stop("The number of elements in group does not match with the number of columns in ymat.")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    print("Conducting the two-sided test as alternative has not been specified...")
  }
  if (alternative=="two.sided"){ side=2 }
  if (alternative=="greater"){ side=1 }
  if (alternative=="less"){ side=-1 }
  if (is.null(seed)){ seed=sample(1e6,1) }
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    
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
        fit=MassiveDiffC(ymatList[[i]], group, nperm, seed, side)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=MassiveDiffC(ymatList[[i]], group, nperm, seed, side)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    return(out)
  } else{
    out=MassiveDiffC(ymat, group, nperm, seed, side)
    out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=TRUE)))/(1+nperm)
    out$seed=seed
    out$alternative=alternative
  }
  return(out)    
}

