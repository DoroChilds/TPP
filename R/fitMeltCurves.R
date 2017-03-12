fitMeltCurves <- function(xMat, yDF, colPrefix, startPars, maxAttempts, expNames, 
                          protID, verbose){
  ## Invoke melting curve fitting for each protein.
  
  if (verbose) message("Fitting melting curves for protein: ", protID)
  
  ## Obtain data for model fit:
  fcCols <- colnames(yDF)[grepl(colPrefix, colnames(yDF))]
  fcMat  <- matrix(nrow=length(expNames), ncol=length(fcCols),
                   dimnames=list(expNames, fcCols))
  fcMat[yDF$expName,fcCols] <- as.matrix(yDF[,fcCols])
  
  listModels <- setNames(vector("list", length(expNames)), expNames)
  
  curveParNames <- meltCurveParamNames(returnParNames = TRUE, 
                                       returnPerformanceInfo = TRUE)
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
      mTmp <- fitSigmoidTR(xVec = xTmp, 
                           yVec = yTmp, 
                           startPars = startPars, 
                           maxAttempts = maxAttempts,
                           fixT0 = TRUE)
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