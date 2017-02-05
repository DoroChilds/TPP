storeMeltCurveParams <- function(data, params){
  ## Store estimated melting curve parameters in featureData.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  expName = protID <- NULL
  
  expNames <- names(data)
  for (en in expNames){
    idsInDataset  <- featureNames(data[[en]])
    
    meltCurveParNames <- c(meltCurveParamNames(returnParNames = TRUE, 
                                               returnPerformanceInfo = TRUE), 
                           "plot")
    
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
