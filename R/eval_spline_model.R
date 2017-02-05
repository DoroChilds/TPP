eval_spline_model <- function(lmFitObj){
  ## Return properties of a given lm object (i.e. from fitting natural smoothing 
  ## splines) like RSS, complexity, data points, sigma, corrected Akaike criterion (AICC),
  ## log-likelihood
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  ranef = fixef <- NULL
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("lmFitObj"))
  
  if (class(lmFitObj) %in% c("lm", "lmerMod")){
    
    ## Residual sum of squares
    rss <- sum(residuals(lmFitObj)^2)
    
    ## Sigma
    sig <- sigma(lmFitObj)
    
    ## corrected Akaike criterion
    aicc <- AICc(lmFitObj)
    
    ## Log-likelihood
    logL <- as.numeric(logLik(lmFitObj))
    
    ## Model matrix for determining number of coefficients and observations
    X <- model.matrix(lmFitObj)
    
    ## Number of observations:
    nObs <- nrow(X)
    
    # Number of coefficients
    if (inherits(lmFitObj, "lmerMod")){
      
      nCoeffsRandom <- length(unlist(ranef(lmFitObj)))
      nCoeffsFixed <- length((fixef(lmFitObj)))
      
      nCoeffs <- nCoeffsRandom + nCoeffsFixed
      
    } else if (inherits(lmFitObj, "lm")){
      
      nCoeffs <- ncol(X)
    }
    
  } else {
    
    rss = nCoeffs = nObs = sig = aicc = logL <- NA
    
  }
  
  fitStats <- data.frame(rss = rss, 
                         nCoeffs = nCoeffs, 
                         nObs = nObs,
                         sigma = sig, 
                         aicc = aicc,
                         loglik = logL)
  
  return(fitStats)
}


