fitMeltCurvesToProtein <- function(xMat, yDF, startPars, maxAttempts, expNames, 
                                   resultPath, plotPathRel, protID, plotTheme, 
                                   grConds, grReps, doPlot, addLegend){
  ## Fit melting curves for each condition and replicate to a single protein
  message("Fitting melting curves for protein: ", protID)
  
  ## Obtain data for model fit:
  fcCols <- colnames(yDF)[grepl("FC", colnames(yDF))]
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
      yTmp <- fcMat[en,]
      flagSufficientDataForFit <- 1
      ## Fit models for current protein:
      mTmp <- fitSigmoidTR(xVec=xTmp, yVec=yTmp, startPars=startPars, 
                           maxAttempts=maxAttempts)
      listModels[[en]] <- mTmp
      if (class(mTmp) != "try-error")   flagModelConverged <- 1
      else                              flagModelConverged <- 0
    } else {
      flagSufficientDataForFit <- 0
      flagModelConverged <- 0
    }
    
    ## Compute melting curve parameters (melting point, inflection point, slope):
    if (flagModelConverged){
      curveParsTmp <- paramsSigmoid(model=mTmp, xRange=range(xTmp, na.rm=TRUE), 
                                    y=yTmp)
    } else{
      curveParsTmp <- rep(NA_real_, (length(curveParNames)-2))
    }
    curveParsWholeProt[en,] <- c(curveParsTmp, flagModelConverged, 
                                 flagSufficientDataForFit)
  }
  
  ## Plot result:
  if (doPlot){
    plotMeltingCurve(modelList=listModels, xMat=xMat, fcMat=fcMat,
                                curvePars=curveParsWholeProt, protID=protID, 
                                filename=file.path(resultPath, plotPathRel), 
                                plotTheme=plotTheme, expConditions=grConds, 
                                expReplicates=grReps, addLegend=addLegend)
  } else {
    plotPathRel <- NA_character_
  }
  
  ## Initialize output:
  curveParsWholeProt$protID <- protID
  curveParsWholeProt$plot   <- plotPathRel
  curveParsWholeProt$expName <- expNames  
  curveParsWholeProt$condition <- grConds
  curveParsWholeProt$replicate <- grReps
  return(curveParsWholeProt)
}