fctSigmoidCCR <- function(){
  ## Return sigmoidal function for curve fit in TPP-CCR experiments
  fctStr <- "1 / (1 + exp((infl - x) * hill))"
  return(fctStr)
}
