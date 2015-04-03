importCheckExperimentNames <- function(expNames, expectedLength){
  ## Assign generic experiment names if not provided by the user
  isGroups <- FALSE
  if (!is.null(expNames)){
    if (length(expNames)==expectedLength){
      isGroups <- TRUE
    }
  }
  if (isGroups == FALSE){
    expNames <- paste("Experiment_", 1:expectedLength, sep="")
    flagDefaultNames <- TRUE
  } else{
    expNames <- as.character(expNames)
    flagDefaultNames <- FALSE    
  }
  out <- list(expNames        = expNames,
              genericExpNames = flagDefaultNames)
  return(out)
}