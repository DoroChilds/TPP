fitMeltCurves <- function(xMat, yDF, colPrefix, startPars, maxAttempts, expNames, 
                          protID, verbose){
  ## Fit melting curves for experiment to a single protein
  if (verbose) message("Fitting melting curves for protein: ", protID)
  
  ## Obtain data for model fit:
  fcCols <- colnames(yDF)[grepl(colPrefix, colnames(yDF))]
  fcMat  <- matrix(nrow=length(expNames), ncol=length(fcCols),
                   dimnames=list(expNames, fcCols))
  fcMat[yDF$expName,fcCols] <- as.matrix(yDF[,fcCols])
  
  listModels <- setNames(vector("list", length(expNames)), expNames)
  
  curveParNames <- meltCurveParamNames()
  curveParsWholeProt<- data.frame(matrix(nrow=length(expNames), 
                                         ncol=length(curveParNames), 
                                         dimnames=list(expNames,curveParNames)))
  
  for (en in expNames){
    xTmp <- xMat[en,]
    yTmp <- fcMat[en,]
    ## Check if number of non-missing values is sufficient
    ## (NLS can only handle data with at least three non-missing values)
    if (sum(!is.na(yTmp))>2){   
      flagSufficientDataForFit <- 1
      ## Fit models for current protein:
      mTmp <- fitSigmoidTR(xVec=xTmp, yVec=yTmp, startPars=startPars, 
                           maxAttempts=maxAttempts)
      listModels[[en]] <- mTmp
      if(class(mTmp) != "try-error"){
        flagModelConverged <- 1
      } else {
        flagModelConverged <- 0
      }
    } else {
      flagSufficientDataForFit <- 0
      flagModelConverged <- 0
    }
    
    curveParsTmp <- rep(NA_real_, (length(curveParNames)-2))
    
    ## Compute melting curve parameters (melting point, inflection point, slope):
    if (flagModelConverged){
      curveParsTmp <- paramsSigmoid(model=mTmp, xRange=range(xTmp, na.rm=TRUE), 
                                    y=yTmp)
      
    }
    
    curveParsWholeProt[en,] <- c(curveParsTmp, flagModelConverged, 
                                 flagSufficientDataForFit)
  } 
  return(list(curveParsWholeProt, fcMat, listModels))
}

fitDRCurve <- function(protID, expName, dose, response, cpd_effect, slBds, verbose){
  ## 1. Preparation:
  flagConv = flagOutsideConcRange = pEC50qualCheckCol <- FALSE
  pec50 = pec50_final = slope = r2 <- NA
  concBds <- guessDRparamBoundaries(dose)
  
  ## 2. Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- sum(!is.na(response))
  flagDat     <- validValues >= 3
  
  ## 3. Perform fit for current protein:
  if (flagDat){
    
    ## Give initial guess on pEC50 & Hill slope
    slBds_tmp <- slBds
    if(cpd_effect=="destabilized") {
      slBds_tmp <- sort(-1*slBds_tmp)
    }    
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
    lbnd <- dose[order(dose)][2]
    ubnd <- dose[order(dose)][length(dose)]
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
  
  outDF <- data.frame("protID"=protID, "pEC50"=pec50_final, "slope"=slope, 
                      "R_sq"=r2, 
                      "pEC50_outside_conc_range"=flagOutsideConcRange, 
                      "pEC50_quality_check" = pEC50qualCheckCol,
                      "model_converged"=flagConv, 
                      "sufficient_data_for_fit"=flagDat,
                      "expName"=expName, stringsAsFactors=FALSE)
  
  return(outDF)
}
