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

buildNNmatrixDist_radius=function(distMat, max.radius=20){
  p=nrow(distMat)
  dist_range = sort(c(ceiling(min(distMat, na.rm=T)),ceiling(min(c(max(distMat,na.rm=T),max.radius),na.rm=T))))
  dist_range_sequence = seq(dist_range[1], dist_range[2], by=1) # different radii to consider
  
  out=foreach(i=1:p,.combine="rbind")%do%{
    dt=unname(distMat[i,])
    ind = which(dt<= dist_range[2]) 
    cbind(i,ind, dt[ind])
  }
  
  NNmatrix=foreach(r=1:length(dist_range_sequence), .combine="rbind")%do%{
    dist=dist_range_sequence[r]
    ind2=which(out[,3]<=dist)
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
  nlocations=do.call("c",lapply(lst, function(x){x$nlocations}))
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
    pvalue=(1+sum(c(permMin)<min(Tstat,na.rm=T)))/(1+nperm)
  } else if (alternative=="greater"){
    threshold=quantile(permMax,1-alpha)
    pvalue=(1+sum(c(permMax)>max(Tstat,na.rm=T)))/(1+nperm)
  } else {
    perm=pmax(abs(permMin),abs(permMax))
    threshold=quantile(perm,1-alpha)
    pvalue=(1+sum(c(perm)>max(abs(Tstat),na.rm=T)))/(1+nperm)
  }

  return(list(
    threshold=threshold,
    Tstat=Tstat,
    permMin=permMin,
    permMax=permMax,
    pvalue=pvalue,
    seed=seed,
    nperm=nperm,
    nlocations=nlocations,
    alternative=alternative
  ))
}


process=function(fit, threshold=NULL){
  if (is.null(threshold)){ threshold=fit$threshold }
  alternative=fit$alternative
  n.locations=fit$nlocations
  cl1=cl2=NULL
  if (alternative=="two.sided"){
    Tmax=apply(matrix(fit$Tstat,n.locations),1,max)
    Tmin=apply(matrix(fit$Tstat,n.locations),1,min)
    cl1=which(Tmax> threshold)
    cl2=which(Tmin< -threshold)
    inter=intersect(cl1,cl2)
    n.inter=length(inter)
    if (n.inter>0){
      for (j in 1:n.inter){
        if (abs(Tmax[inter[j]])>abs(Tmax[inter[j]])){
          cl2=setdiff(cl2, inter[j])
        } else{
          cl1=setdiff(cl1, inter[j])
        }
      }
    }
  } else if (alternative=="greater"){
    cl1=which(apply(matrix(fit$Tstat,n.locations),1,max)> threshold)
  } else if (alternative=="less"){
    cl2=which(apply(matrix(fit$Tstat,n.locations),1,min)< threshold)
  }
  
  return(list(indices.greater=cl1,indices.less=cl2))
}
