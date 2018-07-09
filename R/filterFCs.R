filterFCs <- function(data, pos, lb, ub){
  ## Filter all data sets by specified criteria on fold change columns for 
  ## normalization set construction
  ## Initialize boolean vector indicating valid rows:
  rValid <- rep(TRUE, nrow(data))

  ## Filter proteins according to lower and upper bounds:
  if (max(pos) > ncol(data)){
    msg1 <- "Error during fold change normalization:"
    msg2 <- paste("Given filter criteria require", max(pos), 
                  "fold changes, but only", ncol(data), 
                  "fold change columns are available.")
    msg3 <- "Please adjust the normalization requirements" 
    msg4 <- "(see '?tpptrDefaultNormReqs' for details)."
    stop(paste(msg1, msg2, msg3, msg4))
  }
  for (i in 1:length(pos)){
    rValTmp <- rep(TRUE, nrow(data))

    posTmp <- pos[i]
    lbTmp  <- lb[i]
    ubTmp  <- ub[i]

    x      <- Biobase::exprs(data)[,posTmp]
    lbValid <- x >= lbTmp
    ubValid <- x <= ubTmp

    rValTmp[!lbValid]       <- FALSE
    rValTmp[!ubValid]       <- FALSE
    rValTmp[is.na(lbValid)] <- FALSE
    rValTmp[is.na(ubValid)] <- FALSE
    rValid[rValTmp == FALSE] <- FALSE

    message("  Column ", posTmp, " between ", lbTmp, " and ", ubTmp, " -> ", sum(rValTmp), " out of ", length(rValTmp), " proteins passed")
  }

  fcFiltered = data[rValid,]
  message(sum(rValid), " out of ", length(rValid), " proteins passed in total.\n")
  return(fcFiltered)
}
