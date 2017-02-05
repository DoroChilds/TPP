exportFct_sortCols <- function(dat, idVar, addCol, intensityStr, fcStr, 
                               normalizedData){
  # find sumionarea and different fold change columns
  allCols <- colnames(dat)
  intensityCols  <- grep(intensityStr, allCols, value = TRUE)
  fcOrig         <- grep(paste("^", fcStr, sep=""), allCols, value = TRUE)
  fcTrans        <- grep(".*transformed", allCols, value = TRUE)
  lowestConc     <- grep(".*lowest_conc", allCols, value = TRUE)
  temperatureCol <- grep("temperature", allCols, value = TRUE) # new
  
  # column oder
  if (normalizedData){
    fcNorm <- grep(".*unmodified", allCols, value = TRUE)
    colOrder <- c(idVar, addCol, temperatureCol, intensityCols, fcOrig, fcNorm, 
                  fcTrans)
  }else{
    colOrder <- c(idVar, addCol, temperatureCol, intensityCols, fcOrig, fcTrans)
  }
  restCol <- setdiff(allCols, c(colOrder, lowestConc))
  
  # rearrange columns of dat
  dat <- dat[, c(colOrder, restCol)]
  
  return(dat)
}
