filterOther <- function(data, cols, lb, ub){
  ## Filter all data sets by specified criteria on metadata columns for 
  ## normalization set construction
  
  ## Check if upper and lower bounds are numerics:
  if (!is.numeric(lb) | !is.numeric(ub)){
    stop("Normalization requirements must have numeric upper and lower bounds.")
  }
  
  ## Check if the length of the vectors for upper and lower bounds corresponds 
  ## to the number of desired normalization variables (specified by 'cols'): 
  if (length(cols) != length(lb) | length(cols) != length(ub)){
    longMsg <- paste("Upper bounds ('ub') and lower bounds ('lb') must be",
                      "vectors of the same length as 'cols'.")
    stop(longMsg)
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
      
      longMsg <- paste0("  Column ", colsTmp, " between ", lbTmp, " and ", 
                        ubTmp, "-> ", sum(rValTmp), " out of ", length(rValTmp), 
                        " proteins passed.\n")
      message(longMsg)
    }
  }
  
  for (nm in nonMatches){
    longMsg <- paste("Desired column '",nm,"' not found in the data set. ", 
                 "Therefore, it connot be used as a filter criterion when ", 
                 "creating normalization set.", 
                 sep = "")
    warning(longMsg)
  }
  
  fcFiltered = data[rValid,]
  
  message(sum(rValid)," out of ",length(rValid)," proteins passed in total.\n")
  
  return(fcFiltered)
}
