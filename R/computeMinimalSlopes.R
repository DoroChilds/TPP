computeMinimalSlopes <- function(slV, slT){
  ## Determine the minimum of the melting curve slopes estimated for Vehicle and 
  ## Treatment group for each protein.
  # Conversion will prevent warnings in 'min':
  slV <- mapvalues(slV, NA, Inf, warn_missing=FALSE) 
  slT <- mapvalues(slT, NA, Inf, warn_missing=FALSE) 
  ## Compute minima:
  minSlopes <- apply(cbind(slV, slT), 1, min)
  ## Convert back to prevent Inf values in final output:
  minSlopes <- mapvalues(minSlopes, Inf, NA, warn_missing=FALSE) 
}