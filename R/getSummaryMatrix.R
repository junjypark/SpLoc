getSummaryMatrix=function(ymat, X=NULL, mask,
                          longitudinal=F, n.visits=NULL, randomslope=T, time.var=NULL, 
                          parallel=F, n.cores=1){
  mask.index=which(mask!=0)
  ymat=ymat[mask.index,]
  p=nrow(ymat); n=ncol(ymat)
  if (!longitudinal){
    if (!is.null(X)){
      n=nrow(X); p=ncol(X)
      Q=diag(n)-tcrossprod(tcrossprod(X,solve(crossprod(X))),X)
      ymat=tcrossprod(ymat,Q)
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(n.cores)
      registerDoParallel(cl)
      out=foreach(i=1:p, .combine="rbind")%dopar%{
        ymat[i,]/sum(ymat[i,]^2,na.rm=T)
      }
      stopCluster(cl)
    } else{
      out=ymat
      for (i in 1:p){
        out[i,]=out[i,]/sum(out[i,]^2,na.rm=T)
      }
    }
    
    return(out)
    
  } else{
    print("Fitting a linear mixed model for every voxel/vextex to compare two groups.")
    if (sum(n.visits)!=n){
      stop("The sum of the number of visits does not match with the number of columns of ymat.")
    }

    n.subj=length(n.visits)
    
    timeMat=matrix(0, n.subj, n)
    for (i in 1:n.subj){
      if (i==1){start=1} else{ start=sum(n.visits[1:(i-1)])+1}
      end=sum(n.visits[1:i])
      timeMat[i,start:end]=time[start:end]
    }
    timeMat=Matrix(timeMat, sparse = T)
    
    Subject=rep(paste0("Subj",1:n.subj),n.visits)
    time=X[,time.var]
    lmerctrl=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5))
    if (isTRUE(parallel)){
      cl=makeCluster(n.cores)
      registerDoParallel(cl)
      
      summaryMat=foreach(j=1:p, .combine="rbind", .packages = "SpLoc")%dopar%{
        if (randomslope){
          fit=lmer(ymat[j,] ~ -1+ X+(1+time|Subject), control=lmerctrl)
        } else{
          fit=lmer(ymat[j,] ~ -1+ X+(1|Subject), control=lmerctrl)
        }
        sigma2=attr(VarCorr(fit),"sc")^2
        residuals(fit)/sigma2
      }
      stopCluster(cl)
    } else{
      summaryMat=matrix(NA, p, n)
      for (j in 1:p){
        if (randomslope){
          fit=lmer(ymat[j,] ~ -1+ X+(1+time|Subject),control=lmerctrl)
        } else{
          fit=lmer(ymat[j,] ~ -1+ X+(1|Subject), control=lmerctrl)
        }
        sigma2=attr(VarCorr(fit),"sc")^2
        summaryMat[j,]=residuals(fit)/sigma2
      }
    }
    # out=tcrossprod(summaryMat, timeMat)
    return(list(out=as.matrix(out), TimeMat=timeMat))
  }
}
