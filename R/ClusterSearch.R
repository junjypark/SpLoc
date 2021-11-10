ClusterSearch=function(Tstat, threshold, NNmatrix, alternative=NULL, fraction=0){
  if (length(Tstat)!=NNmatrix@Dim[1]){
    stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
  }
  if (fraction<0 || fraction>1){
    stop("fraction should be between 0 and 1.")
  }
  if (is.null(alternative)){
    print("Assuming two-sided alternative. Please change 'alternative' if not.")
  }
  if (alternative=="less"){
    Tstat=-Tstat
    threshold=-threshold
  }
  
  n.threshold=length(threshold)
  if (n.threshold>1){ threshold=sort(threshold, decreasing = T) } 
  
  selectionList=list()
  TstatList=list()
  NNmatrixList=list()
  for (th in 1:n.threshold){
    if (th==1){ sig.ind=which(Tstat>threshold[th]) } 
    else{ sig.ind=which(Tstat>threshold[th] & Tstat<=threshold[th-1]) }
    TstatList[[th]]=Tstat[sig.ind]
    NNmatrixList[[th]]=NNmatrix[sig.ind,,drop=F]
  }
  
  for (th in 1:n.threshold){
    if (th==1){ selectionList[[th]]=NA }
    else{ selectionList[[th]]=selectionList[[th-1]] }
    
    if (length(TstatList[[th]])==0){ bool=F } else { bool=T }
    while(bool){
      sig=which(NNmatrixList[[th]][which(TstatList[[th]]==max(TstatList[[th]], na.rm=T))[1],]!=0)
      selectionList[[th]]=c(selectionList[[th]],sig)
      
      for (th2 in th:n.threshold){
        num.vec=apply(NNmatrixList[[th2]][,sig,drop=F],1, function(x){sum(x>0)})
        den.vec=as.numeric(table(NNmatrixList[[th2]]@i+1))

        out.set=which(num.vec/den.vec>fraction)
        TstatList[[th2]]=TstatList[[th2]][-out.set]
        NNmatrixList[[th2]]=NNmatrixList[[th2]][-out.set,,drop=F]
      }
      
      if (length(TstatList[[th]])==0){bool=F}
    }
    selectionList[[th]]=c(na.omit(unique(selectionList[[th]])))
  }
  
  return(list(selection=selectionList, threshold=threshold))
}

# prune=function(NegNNmatrix, sig.vec){
#   sig.original=sig=which(sig.vec>0)
#   n.sig=length(sig)
#   if (n.sig>1){
#     NNsig=tcrossprod(NegNNmatrix,sig.vec)
#     ind=which(NNsig>0 & NNsig<n.sig)
#     if (length(ind)>0){
#       NegNNmatrix=NegNNmatrix[ind,,drop=F]
#       ind2=which(NNsig[ind]==as.numeric(table(NegNNmatrix@i+1)))
#       if (length(ind2)==1){
#         exclude=sig[which(NegNNmatrix[ind2,]>0)]
#         sig=setdiff(sig,exclude)
#       } else if (length(ind2)>1){
#         NegNNmatrix=NegNNmatrix[ind2,sig,drop=F]
#         exclude=sig[unique(t(NegNNmatrix)@i+1)]
#         sig=setdiff(sig,exclude)
#       }
#     }
#   }
#   
#   return(list(sig=sig, sig.original=sig.original))
# }

# ClusterSearch2=function(Tstat, threshold, NNmatrix, alternative=NULL, tau=0){
#   if (length(Tstat)!=NNmatrix@Dim[1]){
#     stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
#   }
#   if (tau<0 || tau>1){
#     stop("tau should be between 0 and 1. A recommended value is 0.2")
#   }
#   if (is.null(alternative)){
#     print("Assuming two-sided alternative. Please change 'alternative' if not.")
#   }
#   if (alternative=="less"){
#     Tstat=-Tstat
#     threshold=-threshold
#   }
#   
#   min.threshold=min(threshold)
#   n.threshold=length(threshold)
#   if (n.threshold>1){ threshold=sort(threshold, decreasing = T) } 
#   
#   selectionList=list()
#   TstatList=list()
#   NNmatrixList=list()
#   for (th in 1:n.threshold){
#     if (th==1){ sig.ind=which(Tstat>threshold[th]) } 
#     else{ sig.ind=which(Tstat>threshold[th] & Tstat<=threshold[th-1]) }
#     
#     TstatList[[th]]=Tstat[sig.ind] 
#     NNmatrixList[[th]]=NNmatrix[sig.ind,,drop=F]
#   }
#   
#   ind0=which(Tstat<0)
#   ind1=which(Tstat>min.threshold)
#   ind2=unique(summary(NNmatrix[ind1,])$j)
#   # ind2=unique(summary(NNmatrix[ind1,])[,2])
#   out=unique(NNmatrix[,-ind2]@i+1)
#   NNmatrix=NNmatrix[setdiff(ind0, out),]
#   
#   for (th in 1:n.threshold){
#     if (th==1){ selectionList[[th]]=NA }
#     else{ selectionList[[th]]=selectionList[[th-1]] }
#     
#     if (length(TstatList[[th]])==0){ bool=F } else { bool=T }
#     while(bool){
#       i=which(TstatList[[th]]==max(TstatList[[th]], na.rm=T))[1]
#       sig.prune=prune(NNmatrix, NNmatrixList[[th]][i,,drop=F])
#       sig=sig.prune$sig
#       sig.original=sig.prune$sig.original
#       sig.update=unique(c(selectionList[[th]],sig))
#       selectionList[[th]]=sig.update
#       
#       for (th2 in th:n.threshold){
#         num.vec=apply(NNmatrixList[[th2]][,na.omit(sig.original),drop=F],1,function(x){sum(x>0)})
#         den.vec=as.numeric(table(NNmatrixList[[th2]]@i+1))
#         out.set=unique(c(which(num.vec/den.vec>tau),i))
#         
#         if (length(out.set)>0 & length(TstatList[[th2]]>0)){
#           TstatList[[th2]]=TstatList[[th2]][-out.set]
#           NNmatrixList[[th2]]=NNmatrixList[[th2]][-out.set,,drop=F]
#         }
#       }
#       
#       if (length(TstatList[[th]])==0){bool=F}
#     }
#     selectionList[[th]]=c(na.omit(unique(selectionList[[th]])))
#   }
#   
#   return(list(selection=selectionList, threshold=threshold))
# }

