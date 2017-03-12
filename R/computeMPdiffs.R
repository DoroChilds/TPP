computeMPdiffs <- function(xV, xT){
  ## Determine the difference in melting points estimated for Vehicle and 
  ## Treatment group for each protein.
  mpDiffs <- xT - xV
  return(mpDiffs)
}