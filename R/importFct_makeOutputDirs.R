importFct_makeOutputDirs <- function(outDir, fNames){
  if (is.null(outDir)){
    if(!is.null(fNames)){
      outDir <- dirname(fNames[1])
      outDir <- file.path(outDir, "TPP_results")
      doWrite <- TRUE
    } else {
      m<-"No output directory specified. No result files or plots will be produced."
      message(m)
      doWrite <- FALSE
    }
  } else {
    doWrite <- TRUE
  } 
  
  if (doWrite){
    # Ensure that full path is obtained. (If the path hints at a symbolic link, 
    # there could otherwise be problems with the links embedded in the Excel output):
    if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
    outDir <- normalizePath(outDir, winslash = "/") # Will create warning if directory does not exist yet, therefore we create it first.
    
    message("Results will be written to ", outDir, "\n\n")
    
    ## Create output directory and include a subfolder for data objects created during 
    ## package excecution:
    pathDataObj <- file.path(outDir, "dataObj")
    if (!file.exists(pathDataObj)) dir.create(pathDataObj, recursive=TRUE)
    
  } else {
    pathDataObj <- NULL
  }
  outList <- list(doWrite=doWrite, pathDataObj=pathDataObj, outDir=outDir)
}