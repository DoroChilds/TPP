importCheckTemperatures <- function(temp){
  ## Convert to matrix. Assign rownames manually to ensure consistent output
  ## regardless of 'temp' being native data.frame or tibble.
  tempMatrix <- as.matrix(temp)
  rownames(tempMatrix) <- NULL
  
  ## Make sure that temperature matrix has non-missing values:
  naRows <- apply(is.na(tempMatrix), 1, all)
  
  if(any(naRows)){
    stop("Row(s) ", paste(which(naRows), collapse=", "), 
         " in the configuration table contain only missing temperature values.")
  }
  
  return(tempMatrix)
}
