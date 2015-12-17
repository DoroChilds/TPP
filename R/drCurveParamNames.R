drCurveParamNames <- function(names=TRUE, info=TRUE){
  ##  Assistant function that returns the column names of the dose response 
  ## curve parameters in the internal datasets.
  ## This function is intended to assure consistency when accessing, 
  ## manipulating, or storing melting curve paramter columns in the package's 
  ## data objects.
  out <- c()
  if (names) {
    out <-c(out, "pEC50", "slope", "R_sq")
  }
  if (info) {
    out <- c(out, "pEC50_outside_conc_range", "model_converged", 
             "pEC50_quality_check", "sufficient_data_for_fit")
  }
  return(out)
}