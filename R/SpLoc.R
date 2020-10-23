SpLoc=function(NNmatrix, ymat, nperm=5000, alpha=0.05, seed=NULL){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }

  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }

  if (is.null(seed)){
    print("Specifying a seed value is recommended combining multiple NN matrix.")
  } else{
    set.seed(seed)
  }

  out=SpLocC(NNmatrix, ymat, nperm, alpha)
  out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat)))/(1+nperm)

  return(out)
}

