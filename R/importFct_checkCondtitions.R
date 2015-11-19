importFct_checkConditions <- function(condInfo, expectedLength){
  ## Assign generic condition information to each experiment, if not specified 
  ## by the user.
  flagGenerateConds <- FALSE
  if (is.null(condInfo)){
    message("No information about experimental conditions given.")
    flagGenerateConds <- TRUE
  } else{
    condInfo <- as.character(condInfo)
    condLevels <- unique(condInfo)
    if (!identical(sort(condLevels), c("Treatment", "Vehicle"))){
      message("Condition column contains the values '", paste(sort(condLevels), collapse="', '"), "'.")
      message("Need the values 'Treatment' and 'Vehicle' instead.")
      flagGenerateConds <- TRUE
    }
  }
  if (flagGenerateConds){
    message("Assigning NA to all conditions.
            Reminder: recognition of Vehicle and Treatment groups during pairwise 
            comparisons is only possible when they are specified in the config table.\n")
    condInfo <- rep(NA_character_, expectedLength)
  }
  return(as.character(condInfo))
}
