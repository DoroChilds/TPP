filterOther <- function(data, cols, lb, ub){
  ## Filter all data sets by specified criteria on metadata columns for 
  ## normalization set construction
  
  ## Check if upper and lower bounds are numerics:
  if (!is.numeric(lb) | !is.numeric(ub)){
    stop("Normalization requirements must have numeric upper and lower bounds.")
  }
  
  ## Initialize boolean vector indicating valid rows:
  rValid <- rep(TRUE, nrow(data))
  
  ## Retrieve potential columns for filtering:
  datCols <- featureData(data)@data
  
  ## Filter proteins according to lower and upper bounds:
  matches <- intersect(cols, names(datCols))
  nonMatches <- setdiff(cols, names(datCols))
  if (length(matches) > 0){
    for (i in 1:length(matches)){
      colsTmp <- cols[i]
      lbTmp   <- lb[i]
      ubTmp   <- ub[i]
      
      x       <- datCols[,colsTmp]
      passesLowerBound <- x >= lbTmp
      passesUpperBound <- x <= ubTmp
      
      passesLowerBound[is.na(passesLowerBound)] <- FALSE
      passesUpperBound[is.na(passesUpperBound)] <- FALSE
      
      rValTmp <- passesLowerBound & passesUpperBound
      rValid[!rValTmp] <- FALSE
      
      message("  Column ", colsTmp, " between ", lbTmp, " and ", ubTmp, "-> ", sum(rValTmp), " out of ", length(rValTmp), " proteins passed.\n")
    }
  }
  for (nm in nonMatches){
    msg <- paste("Desired column '",nm,"' not found in the data set. ", 
                 "Therefore, it connot be used as a filter criterion when ", 
                 "creating normalization set.", 
                 sep = "")
    message(msg)
  }
  
  fcFiltered = data[rValid,]
  message(sum(rValid), " out of ", length(rValid), " proteins passed in total.\n")
  return(fcFiltered)
}