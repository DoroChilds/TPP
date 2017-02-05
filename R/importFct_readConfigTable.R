importFct_readConfigTable <- function(cfg){
  ## Check whether obligatory experiment information is provided via 
  ## 1) data frame, or 
  ## 2) spreadsheet file.
  ##
  ## If 2): import from spreadsheet and return table as data frame.
  
  if (is.character(cfg)){
    if (file.exists(cfg)){
      ## Determine file extension and import in appropriate format:
      strChunks <- strsplit(cfg, "\\.")[[1]]
      fileExtension <- strChunks[length(strChunks)]
      if (fileExtension=="txt") {
        ## Import table from tab-delimited file
        tab <- read.table(file=cfg, header=TRUE, check.names=FALSE, 
                          stringsAsFactors=FALSE, sep="\t")
      } else if (fileExtension=="csv"){
        tab <- read.table(file=cfg, header=TRUE, check.names=FALSE, 
                          stringsAsFactors=FALSE, sep=",")
      } else if (fileExtension=="xlsx") {
        ## Import table in Excel format
        tab <- read.xlsx(cfg)
      } else {
        stop("Error during data import: ", cfg, " does not belong to a valid configuration file.")
      }
    } else {
      stop("Error during data import: ", cfg, " does not belong to a valid configuration file.")
    }
    cfg <- tab
  }
  return(cfg)
}


