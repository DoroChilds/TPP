importCheckReplicates <- function(repInfo, expectedLength){
  ## Assign generic replicate information to each experiment, if not specified by the user.
  if (!is.null(repInfo)){
    if (length(repInfo)==expectedLength){
      replValid <- TRUE
    } else {
      message("Vector of experimental replicates has length ", length(repInfo),", which does not match number of experiments (",expectedLength,").")
      replValid <- FALSE
    }
  } else {
    message("No information about experimental replicates given.")
    replValid <- FALSE
  }
  
  if (replValid == TRUE){
    return(as.character(repInfo)) # Convert to character vector so that it replicate numbers can be stored in expressionSet annotation data
  } else {
    message("Assigning value 1 to all replicates.\nReminder: pairwise comparisons between different replicates are currently only possible if exactly two replicates are specified.\n")
    defaultRepl <- rep("1", expectedLength)
    return(defaultRepl)    
  }
}