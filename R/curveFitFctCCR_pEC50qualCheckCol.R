curveFitFctCCR_pEC50qualCheckCol <- function(x, xmin, xmax){
  xInfo <- x
  xInfo[x < xmin] <- "< xmin"
  xInfo[x > xmax] <- NA
  return(xInfo)
}