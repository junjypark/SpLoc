combine=function(lst, alpha=0.05){
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