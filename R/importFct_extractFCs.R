importFct_extractFCs <- function(datDF, colsFC, type){
  ## Extract fold change columns from the imported data tables and return them 
  ## in a matrix.
  
  datFC           <- subset(datDF, select=colsFC)
  datFCChar        <- colwise(as.character)(datFC)
  datFCNum        <- colwise(as.numeric)(datFCChar)
  fcMat           <- as.matrix(datFCNum)
  rownames(fcMat) <- rownames(datFC)
  if (all(is.na(fcMat))){
    stop("Error while loading data: Fold change columns contain only missing values after data import.")
  }
  return(fcMat)
  
}