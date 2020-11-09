getSummaryMatrix=function(ymat, X=NULL, mask,
                          subj.id=NULL, n.visits=NULL, time.var=NULL, random.var=NULL,
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
    Subject=rep(paste0("Subj",1:n.subj),n.visits)
    time=X[,time.var]
    residMat=matrix(NA, p, n.subj)
    for (j in 1:p){
      fit=lmer(ymat[,j] ~ -1+ X+(X[,random.var]|Subject), 
               control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-5)))
      sigma2=attr(VarCorr(fit),"sc")^2
      summaryMat[j,]=time*residuals(fit)/sigma2
    }
    return (summaryMat)
  }
}

