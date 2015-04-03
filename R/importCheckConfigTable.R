importCheckConfigTable <- function(infoTable){
  ## Check whether obligatory experiment information is provided via data 
  #' frame or spreadsheet file.
  if (is.character(infoTable)){
    if (file.exists(infoTable)){
      ## Determine file extension and import in appropriate format:
      strChunks <- strsplit(infoTable, "\\.")[[1]]
      fileExtension <- strChunks[length(strChunks)]
      if (fileExtension=="txt"|fileExtension=="csv") {
        tab <- read.delim(infoTable, check.names=FALSE) ## Import table from tab-delimited file
      } else if (fileExtension=="xlsx") {
        tab <- read.xlsx(infoTable)  ## Import table in Excel format
      } else {
        stop("Error during data import: the name", infoTable, " does not belong to a valid configuration file.")
      }
    } else {
      stop("Error during data import: the name", infoTable, " does not belong to a valid configuration file.")
    }
    infoTable <- tab
  }
  return(infoTable)
}