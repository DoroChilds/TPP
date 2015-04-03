importCheckLabels <- function(labelValues, temperatureMat){
  ## Check label argument for the import function
  if (is.character(labelValues) | is.factor(labelValues)){
    if (length(labelValues) != ncol(temperatureMat)){
      stop("Number of labels in must correspond to number of temperature values per group.")
    }
  } else{
    stop("labels must be a vector of character or factorial variables.")
  }
  labelValues  <- as.character(labelValues)
  return(labelValues)
}