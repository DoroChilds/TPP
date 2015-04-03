storeMeltCurveParams <- function(data, params){
  ## Store estimated melting curve parameters in featureData.
  expNames <- names(data)
  for (en in expNames){
    idsInDataset  <- featureNames(data[[en]])
    meltCurveParNames <- c(meltCurveParamNames(), "plot")
    parsTmp <- subset(params, expName==en & protID %in% idsInDataset, 
                      select=c("protID", meltCurveParNames))
    fillData <- data.frame("protID"=idsInDataset)
    newData <- join(fillData, parsTmp, by="protID")
    pData(featureData(data[[en]]))[,meltCurveParNames] <- newData[,meltCurveParNames]
  }
  return(data)
}