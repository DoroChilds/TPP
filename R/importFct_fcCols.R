importFct_fcCols <- function(datDF, fcPrefix, labelSuffix){
  # Identify fold change columns in the imported data frame.
  # Fold change columns are composed of the string in 'fcPrefix' (prefix of 
  # every fold change column) and the strings contained in 'labelSuffix' (the 
  # TMT labels).
  colsFC    <- paste(fcPrefix, labelSuffix, sep="")
  fcColPos <- match(colsFC, colnames(datDF))
  if (any(is.na(fcColPos))){
    errormsg <- paste("Could not find column(s) '", 
                      paste(colsFC[is.na(fcColPos)], collapse="', '"), 
                      "'. Please specificy the correct fold change column prefix and isotope labels.", 
                      sep="")
    stop(errormsg)
  } else {
    return(colsFC)
  }
}
