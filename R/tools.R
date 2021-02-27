combine=function(lst, alpha=0.05){
  n=length(lst)
  seed=do.call("c",lapply(lst, function(x){x$seed}))
  nperm=do.call("c",lapply(lst, function(x){x$nperm}))
  alternative=do.call("c",lapply(lst, function(x){x$alternative}))
  
  if (length(unique(seed))>1){stop("Use the same seed for every SpLoc output.")}
  else{seed=seed[1]}
  
  if (length(unique(alternative))>1){stop("Use the same alternative for every SpLoc output.")}
  else{alternative=alternative[1]}
  
  if (length(unique(nperm))>1){stop("Use the same number of permutations for every SpLoc output.")}
  else{nperm=nperm[1]}
  
  Tstat=do.call("c",lapply(lst, function(x){x$Tstat}))
  if (alternative=="less"){
    permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,min)
    threshold=quantile(permMax, alpha)
    pvalue=(1+sum(c(permMax)<max(Tstat,na.rm=T)))/(1+nperm[1])
  } else{
    permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)
    threshold=quantile(permMax, 1-alpha)
    pvalue=(1+sum(c(permMax)>max(Tstat,na.rm=T)))/(1+nperm[1])
  }

  return(list(
    threshold=threshold,
    Tstat=Tstat,
    permMax=permMax,
    pvalue=pvalue,
    seed=seed,
    alternative=alternative,
    nperm=nperm
  ))
}

process=function(fit, NNmatrix=NULL, thres){
  if (is.null(NNmatrix)){
    p=length(fit$Tstat)
    index=which(fit$Tstat>thres)
    Tstat=fit$Tstat[index]
    if (length(index)==0){NNmatrix.trim=NULL }
    else{ NNmatrix.trim=sparseMatrix(i=1:length(index), j=index, x=1, dims=c(length(index), p)) }
    return(list(Tstat=Tstat,NNmatrix=NNmatrix.trim))
  } else{
    if (length(fit$Tstat)!=nrow(NNmatrix)){ 
      stop("The number of test statistics must correspond to the rows of the NNmatrix.") 
    }
    index=which(fit$Tstat>thres)
    Tstat=fit$Tstat[index]
    NNmatrix.trim=NNmatrix[index,,drop=F]
    return(list(Tstat=Tstat,NNmatrix=NNmatrix.trim))
  }
}

processSpLocfit=function(names.fit, names.NNmatrix, alpha=0.05, 
                         fit.directory=NULL, NNmatrix.directory=NULL,
                         parallel=F, ncores=1){
  if (length(names.fit)!=length(names.NNmatrix)){
    stop("The numbers of elements in names.fit and names.NNmatrix are not the same.")
  }
  n=length(names.fit)
  
  lst.result=list()
  for (i in 1:n){ lst.result[[i]]=readRDS(paste0(fit.directory,names.fit[i])) }
  result.combine=combine(lst.result, alpha=alpha)
  thres=result.combine$threshold
  
  if (isTRUE(parallel)){
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    lst.thresfit=foreach(i=1:n, .packages=("SpLoc"),.noexport = "SpLocC" )%dopar%{
      if (names.NNmatrix[i]==""){
        NN=NULL
      } else{
        NN=readRDS(paste0(NNmatrix.directory,names.NNmatrix[i]))
      }
      result=lst.result[[i]]
      process(result, NN, min(thres))
      }
    stopCluster(cl)
  }else{
    lst.thresfit=foreach(i=1:n, .packages=("SpLoc"),.noexport = "SpLocC" )%do%{
      if (names.NNmatrix[i]==""){
        NN=NULL
      } else{
        NN=readRDS(paste0(NNmatrix.directory,names.NNmatrix[i]))
      }
      result=lst.result[[i]]
      process(result, NN, min(thres))
    }
  }

  NN=do.call("rbind",lapply(lst.thresfit, function(x){x$NNmatrix}))
  Tstat=do.call("c",lapply(lst.thresfit, function(x){x$Tstat}))
  
  return(list(NNmatrix=NN, Tstat=Tstat, threshold=thres))
}