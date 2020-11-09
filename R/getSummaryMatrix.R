getSummaryMatrix=function(ymat, X=NULL, mask,
                          subj.id=NULL, n.visits=NULL, time.var=NULL, randomslope=T,
                          parallel=F, n.cores=1){
  mask.index=which(mask!=0)
  
  if (is.null(subj.id)){
    if (is.null(X)){
      for ( in mask.index){
        ymat[i,]=ymat[i,]/sum(ymat[i,]^2,na.rm=T)
      }
    } else{
      X=cbind(1,X)
      n=nrow(X); p=ncol(X)
      Q=diag(n)-tcrossprod(tcrossprod(X,solve(crossprod(X))),X)
      ymat=tcrossprod(ymat,Q)
    }
    ymat[is.nan(ymat)]=0
    return(ymat)
  } else{
    print("Fitting a linear mixed model for every voxel/vextex to compare two groups.")
    p=nrow(ymat); n=ncol(ymat)
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
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      
      summaryMat=foreach(i=mask, .combine="rbind", .packages = "lme4")%dopar%{
        if (randomslope){
          fit=lmer(ymat[,j] ~ -1+ X+(1+time|Subject), 
                   control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
        } else{
          fit=lmer(ymat[,j] ~ -1+ X+(1|Subject), 
                   control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
        }
        sigma2=attr(VarCorr(fit),"sc")^2
        residuals(fit)/sigma2
      }
      stopCluster(cl)
    } else{
      summaryMat=matrix(NA, p, n.subj)
      for (j in 1:mask){
        if (randomslope){
          fit=lmer(ymat[,j] ~ -1+ X+(1+time|Subject), 
                   control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
        } else{
          fit=lmer(ymat[,j] ~ -1+ X+(1|Subject), 
                   control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
        }
        sigma2=attr(VarCorr(fit),"sc")^2
        summaryMat[j,]=residuals(fit)/sigma2
      }
    }

    out=tcrossprod(summaryMat, timeMat)
    return(out)
  }
}

