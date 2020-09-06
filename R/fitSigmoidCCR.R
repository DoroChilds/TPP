fitSigmoidCCR <- function(xVec, yVec, hill_init, pec50_init, slopeBounds, 
                          concBounds){
  ## Fit dose response curve to a vector of TPP-CCR measurements.
  
  ## Prepare model fit:
  strSigm <- fctSigmoidCCR()
  fitFct <- as.formula(paste("y ~", strSigm))
  
  ## Attempt model fit by numerical optimization with nls:
  lower <- c(slopeBounds[1], concBounds[1])
  upper <- c(slopeBounds[2], concBounds[2])
  startPars <- list(hill=hill_init, infl=pec50_init)
  m <- try(nls(formula=fitFct, algorithm="port", data=list(x=xVec, y=yVec), 
               start=startPars, lower=lower, upper=upper, na.action=na.exclude),
           silent=TRUE)
  
  ## Check if fit was successful and if estimated parameters have sufficient quality:
  retry <- FALSE
  if(class(m) == "try-error") {
    retry <- TRUE
  } else {
    ## If fit was successful extract pEC50 and Hill slope for quality check
    coeffsTmp <-coef(m) 
    hill  <- coeffsTmp["hill"]
    pec50 <- coeffsTmp["infl"]
    if (!(pec50 >= concBounds[1] & pec50 <= concBounds[2] & sign(hill)==sign(hill_init))){
      retry <- TRUE
    }
  }
  
  ## If fit was not successful, or did not yield satisfactory curve parameters, 
  ## repeat by 'naive' grid search algorithm with nls2:
  if (retry==TRUE){
    startNLS2 <- list(hill=slopeBounds, infl=concBounds)
    # capture output because try with silent option does not work for nls2. Reason: nls2 calls try(nls, ...) without silent option internally.
    cc <- capture.output(type="message",
                         m <- try(nls2(formula=fitFct, algorithm="grid-search", 
                                       data  = list(x=xVec, y=yVec), 
                                       start = startNLS2, na.action=na.exclude)))
  }
  
  return(m)
}
