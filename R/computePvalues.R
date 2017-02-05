computePvalues <- function(minSlopes, mpDiffs, binWidth, type, pAdj){
  ## Compute p-values of melt point differences as described by Cox et al.(2008)
  
  ## Bin melting points according to slope values:
  iValid    <- !is.na(minSlopes)
  minSlopes <- minSlopes[iValid]
  mpDiffs   <- mpDiffs[iValid]
  
  bins <- assignBins(x=minSlopes, w=binWidth)
  
  ## Compute p-values for each bin
  pVals <- rep(NA_real_, length(bins))
  for (b in unique(bins)){
    iBin <- which(bins==b)
    pVals[iBin] <- pValsCurrentBin(mpDiffs[iBin], type=type)
  }
  
  ## Perform Benjamini-Hochberg correction (over all bins)
  pVals <- p.adjust(pVals, pAdj)
  
  ## Output vector
  pOut <- rep(NA_real_, length(minSlopes))
  pOut[iValid] <- pVals
  return(pOut)
}

assignBins <- function(x, w, collapseSmallest = TRUE){
  mpNum       <- length(x)
  binWidthRel <- w/mpNum
  binProp <- sort(unique(c(seq(1, 0, by=-binWidthRel), 0)))
  bounds  <- quantile(x, binProp, na.rm=TRUE)
  bins    <- .bincode(x, bounds, include.lowest=TRUE, right=TRUE)
  
  ## If bin with lowest values is smaller than the others, include it into the
  ## succeeding bin:
  if (collapseSmallest){
    if (sum(bins==1) < w) bins[bins==1] <- 2
  }
  
  return(bins)
}
