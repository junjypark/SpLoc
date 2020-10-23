SpLoc=function(NNmatrix, ymat, nperm=5000, alpha=0.05, seed=NULL){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }

  if (is.null(seed)){
    print("Specifying a seed value is recommended combining multiple NN matrix.")
  }
  
  out=SpLocC(NNmatrix, ymat, nperm, alpha, seed)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat)))/(1+nperm)

  return(out)
}

