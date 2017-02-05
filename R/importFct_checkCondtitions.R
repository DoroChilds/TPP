importFct_checkConditions <- function(condInfo, expectedLength){
  ## Assign generic condition information to each experiment, if not specified 
  ## by the user.
  flagGenerateConds <- FALSE
  
  if (is.null(condInfo)){
    
    message("No information about experimental conditions given. Assigning NA instead.\nReminder: recognition of Vehicle and Treatment groups during pairwise \ncomparisons is only possible when they are specified in the config table.\n")
    
    condInfo <- rep(NA_character_, expectedLength)
    
  } else{
    
    condInfo <- as.character(condInfo) %>% stringr::str_to_title()
    condLevels <- unique(condInfo)
    
    invalidLevels = setdiff(condLevels, c("Treatment", "Vehicle"))
    if (length(invalidLevels) > 0){
      stop("The entry '", invalidLevels, "' in the condition column is invalid. Only the values 'Treatment' and 'Vehicle' are allowed. Please correct this and start again.")
    }
  }
  
  return(condInfo)
}
