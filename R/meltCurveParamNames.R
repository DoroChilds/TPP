meltCurveParamNames <- function(returnParNames=TRUE, returnPerformanceInfo=TRUE){
  ##  Assistant function that returns the column names of the melting curve 
  ## parameters in the internal datasets.
  ## This function is intended to assure consistency when accessing, 
  ## manipulating, or storing melting curve paramter columns in the package's 
  ## data objects.
  out <- c()
  if (returnParNames) {
    out <-c(out, "a", "b", "meltPoint", "inflPoint", "slope", "plateau", "R_sq")
  }
  if (returnPerformanceInfo) {
    out <- c(out, "model_converged", "sufficient_data_for_fit")
  }
  return(out)
}