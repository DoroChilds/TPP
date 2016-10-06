pValFctPerformSingleComparison <- function(minsl, mpdiff, method, control, 
                                           comparisonName){
  if (sum(!is.na(mpdiff)) > 0){
    if (method=="robustZ"){
      ## Compute p-values as described in the MaxQuant paper.
      binWidth = control$binWidth
      nMinSl <- sum(!is.na(minsl))
      if (binWidth > nMinSl){
        warning("P-value computation for comparison ", comparisonName, ":\nAssigned bin width (",binWidth,") is larger than maximum number of minimal slopes available for binning (",nMinSl,"). Setting it to ", nMinSl, " instead.")
        binWidthNew <- nMinSl
      } else{
        binWidthNew <- binWidth
      }
      pVals<-try(computePvalues(minSlopes=minsl, mpDiffs=mpdiff, binWidth=binWidthNew))
      if (class(pVals) == "try-error"){
        pVals <- rep(NA_real_, length(mpdiff))
      }
      
    } else {
      stop("Currently, we can only compute p-values with method 'robustZ'.")
    }
  } else {
    pVals <- rep(NA_real_, length(mpdiff))
    warning("P-value calculation not possible for any protein in comparison '",comparisonName,"'.\nProbably none of the estimated melting curves passed the quality\nchecks with respect to R^2 and plateau (Vehicle).")
  }
  return(pVals)
}