fitDRCurve <- function(protID, expName, dose, response, cpd_effect, slBds, verbose){
  
  ## Preparation:
  flagConv = flagOutsideConcRange = pEC50qualCheckCol <- FALSE
  pec50 = pec50_final = slope = r2 <- NA
  
  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- sum(!is.na(response))
  flagDat     <- validValues >= 3
  
  ## Perform fit for current protein:
  if (flagDat){
    
    ## Sort values:
    o <- order(dose)
    dose <- dose[o]
    response <- response[o]
    
    ## Compute boundaries for optimization:
    concBds <- guessDRparamBoundaries(dose)
    slBds_tmp <- slBds
    if(cpd_effect=="destabilized") {
      slBds_tmp <- sort(-1*slBds_tmp)
    }
    
    ## Give initial guess on pEC50 & Hill slope
    pec50_init <- guessInitialpEC50(dose, response, concBds)
    hill_init  <- guessInitialDRslope(dose, response, slBds_tmp, cpd_effect)
    
    ## Fit curve by numerical optimization (nls and nls2)
    fit <- fitSigmoidCCR(xVec=dose, yVec=response, 
                         hill_init=hill_init, pec50_init=pec50_init,
                         slopeBounds=slBds_tmp, concBounds=concBds)
    
    if(class(fit) != "try-error") {
      # if fit was successful extract pEC50 and Hill slope and calculate R2
      pec50 <- coef(fit)["infl"]
      slope <- coef(fit)["hill"]
      r2    <- rSquared(fit, response)
      pec50_final <- -1*pec50
      flagConv <- TRUE
    }
    
    ## Check if estimated slope is outside of concentration range
    lbnd <- dose[2] # careful: dose must be sorted
    ubnd <- max(dose)
    flagOutsideConcRange <- pec50 > ubnd | pec50 < lbnd
    
    pEC50qualCheckCol <- curveFitFctCCR_pEC50qualCheckCol(x=pec50_final,
                                                          xmin=-ubnd, 
                                                          xmax=-lbnd)
  }
  
  ## Report/ return results:
  if (verbose){
    message("Dose response curve for protein ", protID, ": ", validValues, 
            " non-NA fold changes.\tModel converged: ", flagConv)    
  }
  
  out <- data.frame("protID" = protID, 
                    "pEC50" = pec50_final, 
                    "slope" = slope, 
                    "R_sq" = r2, 
                    "pEC50_outside_conc_range" = flagOutsideConcRange, 
                    "pEC50_quality_check" = pEC50qualCheckCol,
                    "model_converged" = flagConv, 
                    "sufficient_data_for_fit" = flagDat,
                    "expName" = expName, 
                    stringsAsFactors = FALSE)
  
  return(out)
}
