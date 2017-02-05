fitSigmoidTR <- function(xVec, yVec, startPars, maxAttempts, fixT0){
  ## Fit sigmoidal model to a vector of TPP-TR measurements
  strSigm <- fctSigmoidTR(deriv=0)
  fitFct <- as.formula(paste("y ~", strSigm))
  varyPars <- 0
  attempts <- 0
  
  ## For parameter re-sampling in case of non-convergence:
  repeatLoop <- TRUE
  set.seed(123) # make results reproducible
  
  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{  
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.2, 0.2)))      
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), 
                   na.action=na.exclude, algorithm="port", 
                   lower=c(0.0, 1e-5, 1e-5), upper=c(1.5, 15000, 250)), 
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        repeatLoop <- FALSE
      }
    }
  }
  return(m)
}
