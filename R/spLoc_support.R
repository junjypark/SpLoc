
SpLocC <- function(NN, ymat, nperm, alpha){
  .Call('_SpLocC', PACKAGE = "SpLoc", NN, ymat, nperm, alpha)
}
