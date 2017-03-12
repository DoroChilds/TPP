computeMinimalSlopes <- function(xV, xT){
  ## Determine the minimum of the melting curve slopes estimated for Vehicle and 
  ## Treatment group for each protein.
  minSlopes <- pmin(xV, xT)
  return(minSlopes)
}
