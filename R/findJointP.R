findJointP <- function(data){
  ## Find proteins that commonly occur in all experiments. They resulting 
  ## intersect is stored in the variable jointP. It will be used for 
  ## cross-experiment normalization.
  
  message("Detecting intersect between treatment groups (jointP).")
  grNames <- names(data)

  ## Collect protein IDs detected in each treatment group:
  idList <- sapply(grNames, function(n){featureNames(data[[n]])}, simplify=FALSE)

  ## Determine proteins occurring in each table:
  jointP <- idList[[1]]
  for (n in grNames){
    jointP <- intersect(jointP, idList[[n]])
  }
  message("-> JointP contains ", length(jointP), " proteins.\n")
  return(jointP)
}