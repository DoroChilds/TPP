meltCurveParamNames <- function(returnParNames=TRUE, returnPerformanceInfo=TRUE){
  ##  Assistant function that returns the column names of the melting curve 
  ## parameters in the internal datasets.
  ## This function is intended to assure consistency when accessing, 
  ## manipulating, or storing melting curve paramter columns in the package's 
  ## data objects.
  out <- c()
  if (returnParNames) {
    out <-c(out, "a", "b", "meltPoint", "inflPoint", "slope", "plateau", "R_sq")
  }
  if (returnPerformanceInfo) {
    out <- c(out, "model_converged", "sufficient_data_for_fit")
  }
  return(out)
}

drCurveParamNames <- function(names=TRUE, info=TRUE){
  ##  Assistant function that returns the column names of the dose response 
  ## curve parameters in the internal datasets.
  ## This function is intended to assure consistency when accessing, 
  ## manipulating, or storing melting curve paramter columns in the package's 
  ## data objects.
  out <- c()
  if (names) {
    out <-c(out, "pEC50", "slope", "R_sq")
  }
  if (info) {
    out <- c(out, "pEC50_outside_conc_range", "model_converged", 
             "sufficient_data_for_fit")
  }
  return(out)
}

storeMeltCurveParams <- function(data, params){
  ## Store estimated melting curve parameters in featureData.
  expNames <- names(data)
  for (en in expNames){
    idsInDataset  <- featureNames(data[[en]])
    
    meltCurveParNames <- c(meltCurveParamNames(), "plot")
    
    ciOpt <- getOption("TPPTR_CI")
    if(!is.null(ciOpt)){
      if (ciOpt){
        meltCurveParNames <- c(meltCurveParNames, "CI_meltPointUpper", 
                               "CI_meltPointLower", "CI_meltPoint_delta")
      }
    }
    
    parsTmp <- subset(params, expName==en & protID %in% idsInDataset, 
                      select=c("protID", meltCurveParNames))
    fillData <- data.frame("protID"=idsInDataset)
    newData <- join(fillData, parsTmp, by="protID")
    pData(featureData(data[[en]]))[,meltCurveParNames] <- newData[,meltCurveParNames]
  }
  return(data)
}

storeDRCurveParams <- function(data, params){
  ## Store estimated dose response curve parameters in featureData.
  expNames <- names(data)
  curvePars <- drCurveParamNames(TRUE, TRUE)
  for (en in expNames){
    idsInDataset  <- featureNames(data[[en]])
    parsTmp <- subset(params, expName==en & protID %in% idsInDataset, 
                      select=c("protID", curvePars))
    fillData <- data.frame("protID"=idsInDataset)
    newData <- join(fillData, parsTmp, by="protID")
    pData(featureData(data[[en]]))[,curvePars] <- newData[,curvePars]
  }
  return(data)
}
