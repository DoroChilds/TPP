checkIfLegendPossible <- function(){
  ## Determine operating system to decide whether legends should be added.
  osType <- Sys.info()['sysname']
  if (osType %in% c("Linux", "Windows")){
    addLegend <- TRUE    
  } else addLegend <- FALSE
  return(addLegend)
}
