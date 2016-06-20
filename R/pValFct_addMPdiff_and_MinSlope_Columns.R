pValFct_addMPdiff_and_MinSlope_Columns <- function(expNameV, expNameT, parDF){
  ## Helper function for p-value calculation:
  ## Compute melting point differences and minimal slopes, and attach to data 
  ## frame with fitted melting curve parameters.
  
  mpV <- parDF[, paste("meltPoint", expNameV, sep="_")]
  mpT <- parDF[, paste("meltPoint", expNameT, sep="_")]
  slV <- parDF[, paste("slope", expNameV, sep="_")]
  slT <- parDF[, paste("slope", expNameT, sep="_")]
  mpDiffs <- computeMPdiffs(xV=mpV, xT=mpT)
  minSl   <- computeMinimalSlopes(xV=slV, xT=slT)
  
  ## Store computed values in output table:
  mpDiffDF <- data.frame(mpDiffs, stringsAsFactors=FALSE)
  minSlDF   <- data.frame(minSl, stringsAsFactors=FALSE)
  colnames(mpDiffDF) <- paste("diff_meltP", expNameT, "vs", expNameV, sep="_")
  colnames(minSlDF)  <- paste("min_slope" , expNameT, "vs", expNameV, sep="_")
  
  return(list(mpDiffs=mpDiffDF,
              minSl  = minSlDF))
  
}

computeMinimalSlopes <- function(xV, xT){
  ## Determine the minimum of the melting curve slopes estimated for Vehicle and 
  ## Treatment group for each protein.
  
  minSlopes <- pmin(xV, xT)
  return(minSlopes)
}

computeMPdiffs <- function(xV, xT){
  ## Determine the difference in melting points estimated for Vehicle and 
  ## Treatment group for each protein.
  
  mpDiffs <- xT - xV
  return(mpDiffs)
}