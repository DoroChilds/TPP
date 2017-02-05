removeInvalidPaths <- function(paths){
  # Remove non-existing files:
  # Problem: it is highly inefficient to try the on every element, even if 
  # plotting was not invoked at all (approx. 1 minute for 60,000 file names)
  isValidPath <- file.exists(paths) 
  validPaths  <- ifelse(isValidPath, yes = paths, no = NA)
  
  # Remove potential duplicates:
  isDuplicated <- duplicated(validPaths, imcomparables = NA)
  uniquePaths  <- ifelse(isDuplicated, yes = NA, no = validPaths)
  
  return(uniquePaths)
  
  # createAndCheckPaths <- function(prefixVec, suffix){
  #   
  #   ## Create plot paths:
  #   allPaths <- paste0(prefixVec, suffix)
  #   
  #   ## Check files for existence and remove duplicates:
  #   validPaths <- removeInvalidPaths(allPaths)
  #   
  #   return(validPaths)
  # }
}

