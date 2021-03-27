ClusterSearch=function(Tstat, threshold, NNmatrix, alternative=NULL, fraction=0){
  if (length(Tstat)!=nrow(NNmatrix@Dim[1])){
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

ClusterSearch2=function(Tstat, threshold, NNmatrix, alternative=NULL, tau=0, is.eps=TRUE, eps=NULL){
  if (length(Tstat)!=nrow(NNmatrix@Dim[1])){
    stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
  }
  if (tau<0 || tau>1){
    stop("tau should be between 0 and 1. A recommended value is 0.2")
  }
  if (is.null(alternative)){
    print("Assuming two-sided alternative. Please change 'alternative' if not.")
  }
  if (alternative=="less"){
    Tstat=-Tstat
    threshold=-threshold
  }
  if (is.eps & is.null(eps)){
    if (alternative=="two.sided"){
      eps=1; print("eps set as 1 as a default")      
    } else {
      eps=0; print("eps set as 0 as a default")
    }
  }
  
  n.threshold=length(threshold)
  if (n.threshold>1){ threshold=sort(threshold, decreasing = T) } 
  
  if (is.eps){
    NNmatEps=NNmatrix[which(Tstat<eps),]
    np=length(NNmatEps@p)
    exclude=which(NNmatEps@p[-1]-NNmatEps@p[-np]>0)
    include=setdiff(1:ncol(NNmatrix), exclude) 
  }
  
  selectionList=list()
  TstatList=list()
  NNmatrixList=list()
  for (th in 1:n.threshold){
    if (th==1){ sig.ind=which(Tstat>threshold[th]) } 
    else{ sig.ind=which(Tstat>threshold[th] & Tstat<=threshold[th-1]) }
    
    if (is.eps){
      TstatTemp=Tstat[sig.ind]
      NNtemp=NNmatrix[sig.ind,,drop=F]
      NNtemp[,exclude]=0
      ind=sort(unique(NNtemp@i+1))
      TstatList[[th]]=Tstat[ind]
      NNtemp=NNtemp[ind,]
      NNmatrixList[[th]]=NNtemp
    } else{
      TstatList[[th]]=Tstat[sig.ind] 
      NNmatrixList[[th]]=NNmatrix[sig.ind,,drop=F]
    }
  }
  
  for (th in 1:n.threshold){
    if (th==1){ selectionList[[th]]=NA }
    else{ selectionList[[th]]=selectionList[[th-1]] }
    
    if (length(TstatList[[th]])==0){ bool=F } else { bool=T }
    while(bool){
      i=which(TstatList[[th]]==max(TstatList[[th]], na.rm=T))[1]
      sig=which(NNmatrixList[[th]][i,]!=0)
      if (is.eps){ sig=intersect(sig, include) }
      
      selectionList[[th]]=c(selectionList[[th]],sig)
      print(paste0("search done: ",length(sig)))
      
      for (th2 in th:n.threshold){
        num.vec=apply(NNmatrixList[[th2]][,sig,drop=F],1, function(x){sum(x>0)})
        den.vec=as.numeric(table(NNmatrixList[[th2]]@i+1))
        out.set=unique(c(which(num.vec/den.vec>tau),i))
        
        if (length(out.set)>0){
          TstatList[[th2]]=TstatList[[th2]][-out.set]
          NNmatrixList[[th2]]=NNmatrixList[[th2]][-out.set,,drop=F]
        }
      }
      
      if (length(TstatList[[th]])==0){bool=F}
    }
    selectionList[[th]]=c(na.omit(unique(selectionList[[th]])))
  }
  
  return(list(selection=selectionList, threshold=threshold))
}

# ClusterSearchOld=function(Tstat, threshold, NNmatrix){
#   if (length(Tstat)!=nrow(NNmatrix)){
#     stop("The number of rows in NNmat needs to be the same as the length of Tstat.")
#   }
#   
#   if (length(threshold)>1){ threshold=sort(threshold,decreasing = T) } 
#   
#   selection=list()
#   for (th in 1:length(threshold)){
#     if (th==1){ sig.ind=which(Tstat>threshold[th]) } 
#     else{ sig.ind=which(Tstat>threshold[th] & Tstat<=threshold[th-1]) }
#     
#     Tstat.sub=Tstat[sig.ind]
#     NNmatrix.sub=NNmatrix[sig.ind,,drop=F]
#     
#     if (th==1){
#       sig=NULL
#     } else{
#       if (!is.null(sig)){
#         out.set=which(apply(NNmatrix.sub[,sig,drop=F],1,max)>0)
#         if (length(out.set)>0){
#           Tstat.sub=Tstat.sub[-out.set]
#           NNmatrix.sub=NNmatrix.sub[-out.set,,drop=F]
#         }
#       }
#     }
#     
#     if (length(Tstat.sub)==0){bool=F}else{bool=T}
#     
#     while(bool){
#       ind=which(Tstat.sub==max(Tstat.sub))[1]
#       sig.vertices=which(NNmatrix.sub[ind,]!=0)
#       sig=c(sig, sig.vertices)
#       out.set=which(apply(NNmatrix.sub[,sig.vertices,drop=F],1,max)>0)
#       
#       Tstat.sub=Tstat.sub[-out.set]
#       NNmatrix.sub=NNmatrix.sub[-out.set,,drop=F]
#       if (length(Tstat.sub)==0){bool=F}
#     }
#     
#     if (is.null(sig)){ selection[[th]]=NA }
#     else{ selection[[th]]=sig }
#   }
#   
#   return(list(selection=selection, threshold=threshold))
# }