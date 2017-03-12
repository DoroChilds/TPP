pValsCurrentBin <- function(mpDiffs, type){
  ## Compute p-values for each melting point difference.
  
  ## Determine Tm values outside of 0.16 and 0.86 percentiles for each bin.
  mpPerc <- quantile(mpDiffs, probs = c(0.1587, 0.5, 0.8413), na.rm=TRUE)
  r_1 = mpPerc[1]
  r0 = mpPerc[2]
  r1 = mpPerc[3]
  pVals = apply(as.matrix(mpDiffs), MARGIN=1, FUN=mpSinglePval, r_1=r_1, r0=r0, 
                r1=r1, type=type)
  return(pVals)
}