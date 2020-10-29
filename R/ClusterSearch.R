ClusterSearch=function(Tstat, threshold, NNmatrix){
  if (length(Tstat)!=nrow(NNmatrix)){
    stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
  }
  
  sig.ind=which(Tstat>threshold)
  Tstat=Tstat[sig.ind]
  NNmatrix=NNmatrix[sig.ind,,drop=F]
  
  if (nrow(NNmatrix)==0){bool=F}else{bool=T}
  
  clust=1
  sig=NULL
  while(bool){
    ind=which(Tstat==max(Tstat))[1]
    sig.vertices=which(NNmatrix[ind,]!=0)
    sig=c(sig, sig.vertices)
    out.set=which(apply(NNmatrix[,sig.vertices,drop=F],1,max)>0)
    
    Tstat=Tstat[-out.set]
    NNmatrix=NNmatrix[-out.set,,drop=F]
    
    print(paste0("Cluster ", clust, " with ", length(sig.vertices), " voxel(s) is selected (total=", length(sig), ")." ))
    clust=clust+1
    
    if (length(Tstat)==0){bool=F}
  }
  
  return(sig)
}

