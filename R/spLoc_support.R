
SpLocC <- function(NN, ymat, nperm, alpha){
  .Call('_SpLoc_SpLocC', PACKAGE = "SpLoc", NN, ymat, nperm, alpha)
}
