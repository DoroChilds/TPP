filterOther <- function(data, cols, lb, ub){
  ## Filter all data sets by specified criteria on metadata columns for 
  ## normalization set construction
  ## Initialize boolean vector indicating valid rows:
  rValid <- rep(TRUE, nrow(data))
  
  ## Retrieve potential columns for filtering:
  datCols <- featureData(data)@data
  
  ## Filter proteins according to lower and upper bounds:
  for (i in 1:length(cols)){
    if( cols[i] %in% names(datCols) ){
      rValTmp <- rep(TRUE, nrow(data))
      
      colsTmp <- cols[i]
      lbTmp   <- lb[i]
      ubTmp   <- ub[i]
      
      x       <- datCols[,colsTmp]
      lbValid <- x >= lbTmp
      ubValid <- x <= ubTmp
      
      rValTmp[!lbValid]       <- FALSE
      rValTmp[!ubValid]       <- FALSE
      rValTmp[is.na(lbValid)] <- FALSE
      rValTmp[is.na(ubValid)] <- FALSE
      rValid[rValTmp == FALSE] <- FALSE
      
      message("  Column ", colsTmp, " between ", lbTmp, " and ", ubTmp, "-> ", sum(rValTmp), " out of ", length(rValTmp), " proteins passed")
      
    }else{
      warning(paste( cols[i] , "not found in input columns, can't be used for filtering."))
    }
    
  }
  
  fcFiltered = data[rValid,]
  message(sum(rValid), " out of ", length(rValid), " proteins passed in total.\n")
  return(fcFiltered)
}