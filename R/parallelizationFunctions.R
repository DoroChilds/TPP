checkIfLegendPossible <- function(){
  ## Determine operating system to decide whether legends should be added.
  
  osType <- Sys.info()['sysname']
  if (osType %in% c("Linux", "Windows")){
    addLegend <- TRUE    
  } else addLegend <- FALSE
  return(addLegend)
}

checkCPUs <- function(cpus) {
  ## Determine appropriate number of CPUs for parallelization.
  
  maxCores <- detectCores()
  if (is.numeric(cpus) && (cpus > maxCores)) {
    warning("Selected number of cores (", cpus, ") exceeds those available on your device (", maxCores, "). Using ", maxCores, " cores.")
    cpus <- maxCores
  } else if (identical(cpus, "max")) {
    cpus <- maxCores
  } else if (!is.numeric(cpus)){
    stop(paste("Invalid argument",cpus, "for 'cpus'."))
  }
  return(cpus)
}

stopImplicitCluster <- function(){
  .options <- doParallel:::.options 
  if(exists(".revoDoParCluster", where=.options) && 
       !is.null(.options[['.revoDoParCluster']]))
  {
    stopCluster(.options[['.revoDoParCluster']])
    remove('.revoDoParCluster', envir=.options)
  }
}
