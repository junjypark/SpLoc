ClusterSearch=function(tstat, thres, NNmatrix){
  if (length(tstat)!=nrow(NNmatrix)){
    stop("The number of rows in NNmat needs to be the same as the length of tstat.")
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
    
    print(paste0("Cluster ", clust, " with ", length(sig.vertices), " voxel(s) is selected (total=", length(sig), ")." ))
    clust=clust+1
    
    if (length(tstat)==0){bool=F}
  }
  
  return(sig)
}

