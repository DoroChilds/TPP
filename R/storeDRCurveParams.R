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