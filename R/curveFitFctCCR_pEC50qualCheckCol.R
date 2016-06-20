curveFitFctCCR_pEC50qualCheckCol <- function(x, xmin, xmax){
  ## Check whether the pEC50 values calculated from dose response curves lie
  ## within a predefined region. Annotate values outside this regiom for 
  ## reporting in the result table.
  xInfo <- x
  xInfo[x < xmin] <- "< xmin"
  xInfo[x > xmax] <- NA
  return(xInfo)
}