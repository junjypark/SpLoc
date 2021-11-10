buildNNmatrix3D=function(template, radiusSet=c(0,1,2,3)){

  radiusSet=sort(radiusSet)
  nSet=length(radiusSet)
  maxRadius=radiusSet[nSet]
  ind=which(template!=0)
  dim.template=dim(template)
  index.template=which(template!=0,arr.ind = T)
  index=which(template!=0)
  p=sum(template!=0)
  
  v=NULL; distvec=NULL
  for (i in -maxRadius:maxRadius){
    for (j in -maxRadius:maxRadius){
      for (k in -maxRadius:maxRadius){
        dist=sqrt(i^2+j^2+k^2)
        if (dist<=maxRadius){
          v=rbind(v, c(i,j,k))
          distvec=c(distvec, dist)
        }
      }
    }
  }
  
  nv=nrow(v)
  mat=matrix(NA, nv, p)
  
  for (i in 1:p){
    pos=index.template[i,]
    mat[,i]=array2vec(v+rep(pos, each = nv),dim.template)
  }
  
  ivec=kronecker(1:p, rep(1,nv))
  jvec=c(mat)
  distvec2=rep(distvec, p)
  
  NNmatrix=foreach(r=1:nSet, .combine="rbind")%do%{
    radius=radiusSet[r]
    phi=phiSet[r]
    if (radius==0){
      spDist=sparseMatrix(i=1:length(index),j=1:length(index),x=1, dims=c(length(index),length(index)))
    } else{
      index2=which(jvec>0 & distvec2<=radius)
      if (is.null(kernel)){
        spDist=sparseMatrix(i=ivec[index2],j=jvec[index2],x=1, dims=c(p,prod(dim.template)))
      } else if (kernel=="Gaussian"){
        spDist=sparseMatrix(i=ivec[index2],j=jvec[index2],x=exp(-distvec2[index2]^2/phi), dims=c(p,prod(dim.template)))
      } else if (kernel=="Exponential"){
        spDist=sparseMatrix(i=ivec[index2],j=jvec[index2],x=exp(-distvec2[index2]/phi), dims=c(p,prod(dim.template)))
      }
      spDist=spDist[,index]
    }
    spDist
  }
  
  return(NNmatrix)
}

buildNNmatrixDist=function(distMat, nnSet=c(1,5,10*1:10,50*3:10,100*6:10)){
  
  p=nrow(distMat)
  nnSet=unique(pmin(sort(nnSet),p))
  n.nnSet=length(nnSet)
  nnMax=nnSet[n.nnSet]
  
  out=foreach(i=1:p,.combine="rbind")%do%{
    dt=distMat[i,]
    rk=rank(dt)
    ind=which(rk<=nnMax)
    cbind(i,ind, rk[ind], dt[ind])
  }
  
  NNmatrix=foreach(r=1:n.nnSet, .combine="rbind")%do%{
    nn=nnSet[r]
    phi=phiSet[r]
    ind2=which(out[,3]<=nn)
    out2=out[ind2,]
    sp=sparseMatrix(i=out2[,1], j=out2[,2], x=1, dims=c(p,p))
    sp
  }
  
  return(NNmatrix)
}

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
  permMin=apply(do.call("cbind",lapply(lst, function(x){x$permMin})),1,min)
  permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)

  if (alternative=="less"){
    threshold=quantile(permMin,alpha)
    pvalue=(1+sum(c(permMin)<min(Tstat,na.rm=T)))/(1+nperm[1])
  } else if (alternative=="greater"){
    threshold=quantile(permMax,1-alpha)
    pvalue=(1+sum(c(permMax)>max(Tstat,na.rm=T)))/(1+nperm[1])
  } else {
    perm=pmax(abs(permMin),abs(permMax))
    threshold=quantile(pmax(abs(permMin),abs(permMax)),1-alpha)
    pvalue=(1+sum(c(perm)>max(abs(Tstat),na.rm=T)))/(1+nperm[1])
  }

  return(list(
    threshold=threshold,
    Tstat=Tstat,
    permMin=permMin,
    permMax=permMax,
    pvalue=pvalue,
    seed=seed,
    nperm=nperm,
    alternative=alternative
  ))
}

process=function(fit, threshold=NULL){
  if (is.null(threshold)){ threshold=fit$threshold }
  alternative=fit$alternative
  n.locations=fit$n.locations
  cl1=cl2=NULL
  if (alternative=="two.sided"){
    cl1=which(apply(matrix(fit$Tstat,n.locations),1,max)> thres)
    cl2=which(apply(matrix(fit$Tstat,n.locations),1,min)< -thres)
    inter=intersect(cl1,cl2)
    n.inter=length(inter)
    if (n.inter>0){
      for (j in 1:n.inter){
        if (fit$Tstat[inter[j]]>0){
          cl2=setdiff(cl2, inter[j])
        } else{
          cl1=setdiff(cl1, inter[j])
        }
      }
    }
  } else if (alternative=="greater"){
    cl1=which(apply(matrix(fit$Tstat,n.locations),1,max)> thres)
  } else if (alternative=="less"){
    cl2=which(apply(matrix(fit$Tstat,n.locations),1,min)< thres)
  }
  
  return(list(indices.greater=cl1,indices.less=cl2))
}

# process=function(fit, NNmatrix=NULL, thres){
#   if (is.null(NNmatrix)){
#     p=length(fit$Tstat)
#     index=which(fit$Tstat>thres)
#     Tstat=fit$Tstat[index]
#     if (length(index)==0){NNmatrix.trim=NULL }
#     else{ NNmatrix.trim=sparseMatrix(i=1:length(index), j=index, x=1, dims=c(length(index), p)) }
#     return(list(Tstat=Tstat,NNmatrix=NNmatrix.trim))
#   } else{
#     if (length(fit$Tstat)!=nrow(NNmatrix)){ 
#       stop("The number of test statistics must correspond to the rows of the NNmatrix.") 
#     }
#     index=which(fit$Tstat>thres)
#     Tstat=fit$Tstat[index]
#     NNmatrix.trim=NNmatrix[index,,drop=F]
#     return(list(Tstat=Tstat,NNmatrix=NNmatrix.trim))
#   }
# }
