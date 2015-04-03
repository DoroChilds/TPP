importCheckTemperatures <- function(temp, nData){
  ## Determine matrix of temperatures corresponding to the TMT labels
  isTempMat <- FALSE
  isTempVec <- FALSE
  if (is.matrix(temp) || is.data.frame(temp)){
    if (nrow(temp) == nData) isTempMat <- TRUE
  }
  if (!isTempMat){
    if (is.numeric(temp)) isTempVec <- TRUE
  }
  if (isTempMat){
    tempMatrix <- as.matrix(temp)
  } else if (isTempVec){
    tempMatrix <- matrix(temp, nrow=nData, ncol=length(temp), byrow=TRUE)
  } else {
    stop("temperature must either be a numeric vector, or a matrix with one row per treatment group.")
  }
  return(tempMatrix)
}
