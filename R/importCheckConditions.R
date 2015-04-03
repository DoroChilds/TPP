importCheckConditions <- function(condInfo, expectedLength){
  ## Assign generic condition information to each experiment, if not specified 
  ## by the user.
  if (!is.null(condInfo)){
    if (length(condInfo)==expectedLength){
      condsValid <- TRUE
    } else {
      condsValid <- FALSE
      message("Vector of experimental conditions has length ", length(condInfo),", which does not match number of experiments (",expectedLength,").")
    }
  } else {
    condsValid <- FALSE
    message("No information about experimental conditions given.")
  }
  
  if (condsValid == TRUE){
    return(as.character(condInfo))
  } else {
    message("Assigning NA to all conditions.\nReminder: pairwise comparisons between Vehicle and Treatment are only possible if the conditions are specified correctly.\n")
    defaultConds <- rep(NA_character_, expectedLength)
    return(defaultConds)    
  }
}