importFct_checkExperimentCol <- function(expCol){
  if (is.null(expCol)){
    m <- "Config table needs an 'Experiment' column with unique experiment IDs."
    stop(m, "\n")
  }
  
  oldExpNames <- expCol
  newExpNames <- gsub("([^[:alnum:]])", "_", expCol)
  iChanged <- oldExpNames != newExpNames
  if (any(iChanged)){
    m1 <- "Replaced non-alphanumeric characters in the 'Experiment' column entries:"
    m2 <- paste("'", paste(oldExpNames[iChanged], collapse="', '"),
                "'\nby\n'", 
                paste(newExpNames[iChanged], collapse="', '"), sep="")
    message(m1,"\n",m2, "\n")
  }
  return(newExpNames)
  
}