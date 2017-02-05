importFct_create_pData <- function(labels, labelValues, fcCols, type){
  if (type == "TR"){
    dataTmp <- data.frame("label"=labels, "temperature"=labelValues, 
                          "normCoeff"=NA, stringsAsFactors=FALSE, 
                          row.names=fcCols)
    metaTmp   <- data.frame(labelDescription = c("Isobaric label",
                                                 "Temperature",
                                                 "Applied normalization coefficient"),
                            row.names = colnames(dataTmp))
  } else if (type == "CCR"){
    dataTmp <- data.frame("label"=labels, "concentration"=labelValues, 
                          "normCoeff"=NA, stringsAsFactors=FALSE, 
                          row.names=fcCols)
    metaTmp   <- data.frame(labelDescription = c("Isobaric label",
                                                 "Concentration",
                                                 "Applied normalization coefficient"),
                            row.names = colnames(dataTmp))
  }
  colInfo <- AnnotatedDataFrame(data=dataTmp, varMetadata=metaTmp)
  return(colInfo)
}
