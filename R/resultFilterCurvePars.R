resultFilterCurvePars <- function(r2V, r2T, plV, minR2=0.8, maxPlateau=0.3){
  ## Check melting curve quality before starting p-value computation.
  passedTest <- r2V>minR2 & r2T>minR2 & plV<maxPlateau
  passedTest[is.na(passedTest)] <- FALSE
  return(passedTest)
}