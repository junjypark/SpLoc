ClusterSearch=function(Tstat, threshold, NNmatrix){
  if (length(Tstat)!=nrow(NNmatrix)){
    stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
  }
  
  if (length(threshold)>1){ threshold=sort(threshold,decreasing = T) } 
  
  selection=list()
  for (th in 1:length(threshold)){
    if (th==1){ sig.ind=which(Tstat>threshold[th]) } 
    else{ sig.ind=which(Tstat>threshold[th] & Tstat<=threshold[th-1]) }
    
    Tstat.sub=Tstat[sig.ind]
    NNmatrix.sub=NNmatrix[sig.ind,,drop=F]
    
    if (th==1){
      sig=NULL
    } else{
      if (!is.null(sig)){
        out.set=which(apply(NNmatrix.sub[,sig,drop=F],1,max)>0)
        if (length(out.set)>0){
          Tstat.sub=Tstat.sub[-out.set]
          NNmatrix.sub=NNmatrix.sub[-out.set,,drop=F]
        }
      }
    }
    
    if (length(Tstat.sub)==0){bool=F}else{bool=T}
    
    while(bool){
      ind=which(Tstat.sub==max(Tstat.sub))[1]
      sig.vertices=which(NNmatrix.sub[ind,]!=0)
      sig=c(sig, sig.vertices)
      out.set=which(apply(NNmatrix.sub[,sig.vertices,drop=F],1,max)>0)
      
      Tstat.sub=Tstat.sub[-out.set]
      NNmatrix.sub=NNmatrix.sub[-out.set,,drop=F]
      if (length(Tstat.sub)==0){bool=F}
    }
    
    if (is.null(sig)){ selection[[th]]=NA }
    else{ selection[[th]]=sig }
  }
  
  return(list(selection=selection, threshold=threshold))
}


ClusterSearch2=function(Tstat, threshold, NNmatrix, fraction=0.1){
  
  selection=NULL
  sig.ind=which(Tstat>threshold)  
  if (length(sig.ind)>0){
    Tstat.sub=Tstat[sig.ind]
    NNmatrix.sub=NNmatrix[sig.ind,,drop=F]
    
    if (length(Tstat.sub)==0){bool=F}else{bool=T}
    while(bool){
      sig=which(NNmatrix.sub[which(Tstat.sub==max(Tstat.sub, na.rm=T))[1],]!=0)
      selection=c(selection,sig)
      num.vec=apply(NNmatrix.sub[,sig,drop=F],1, function(x){sum(x>0)})
      
      b=which(NNmatrix.sub>0,arr.ind=T)
      den.vec=as.numeric(table(b[,1]))
      
      out.set=which(num.vec/den.vec>fraction)
      Tstat.sub=Tstat.sub[-out.set]
      NNmatrix.sub=NNmatrix.sub[-out.set,,drop=F]
      if (length(Tstat.sub)==0){bool=F}
    }
    selection=unique(selection)
  }
  
  return(list(selection=selection, threshold=threshold))
}

# 
# Booster=function(fit, NNmatrix, parallel=F, ncores=1){
#   ind=which(fit$Tstat>fit$thres)
#   Tstatsub=fit$Tstat[ind]
#   NNsub=NNmatrix[ind,]
#   nonzero.index=which(NNsub!=0, arr.ind=T)
#   voxels=sort(unique(nonzero.index[,2]))
#   
#   if (parallel){
#     cl=makeCluster(ncores)
#     registerDoParallel(cl)
#     boost=foreach(i=1:length(voxels), .combine="c", .packages="Matrix")%dopar%{
#       max(Tstatsub[which(NNsub[,voxels[i]]!=0)])
#     }
#     stopCluster(cl)
#   } else{
#     boost=foreach(i=1:length(voxels), .combine="c",.packages="Matrix")%do%{
#       max(Tstatsub[which(NNsub[,voxels[i]]!=0)])
#     }
#   }
#   
#   return(list(boost=boost, voxels=voxels))
# }
# 
