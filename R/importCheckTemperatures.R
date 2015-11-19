importCheckTemperatures <- function(temp){
  tempMatrix <- as.matrix(temp)
  ## Make sure that temperature matrix has non-missing values:
  naRows <- apply(is.na(tempMatrix), 1, all)
  
  if(any(naRows)){
    stop("Row(s) ", paste(which(naRows), collapse=", "), 
         " in the configuration table contain only missing temperature values.")
  }
  
  return(tempMatrix)
}
