importFct_readFiles <- function(files, naStrs){
  ## Import a single TPP-TR or -CCR experiment from files
  expNames    <- names(files)
  data        <- vector("list", length(files))
  names(data) <- expNames
  
  for (expName in expNames){
    fTmp <- files[[expName]]
    if (file.exists(fTmp) || url.exists(fTmp)){
      data[[expName]] <- read.delim(fTmp, as.is=TRUE, na.strings=naStrs, quote="")      
    } else {
      stop("File ", fTmp, " could not be found.")
    }
  }
  return(data)
  
  # to do: incorporate the following code:
  # importFct_2Ddataframe <- function(filePath, rowNumber){
  #   if (file.exists(filePath)){
  #     ## Determine file extension and import in appropriate format:
  #     strChunks <- strsplit(filePath, "\\.")[[1]]
  #     fileExtension <- strChunks[length(strChunks)]
  #     if (fileExtension=="txt") {
  #       ## Import table from tab-delimited file
  #       tab <- read.delim(file=filePath, header=TRUE, check.names=FALSE, as.is=TRUE,
  #                         stringsAsFactors=FALSE, sep="\t", blank.lines.skip = FALSE)
  #     } else if (fileExtension=="csv"){
  #       tab <- read.table(file=filePath, header=TRUE, check.names=FALSE, as.is=TRUE,
  #                         stringsAsFactors=FALSE, sep=",", blank.lines.skip = FALSE)
  #     } else if (fileExtension=="xlsx") {
  #       ## Import table in Excel format
  #       tab <- read.xlsx(filePath)
  #     } else {
  #       stop(paste("Error during data import: ", filePath, " in line ", as.character(rowNumber), " 
  #                  of the configTable does not correspond to not an existing file!", sep=""))
  #     }
  #     return(tab)
  #     } else{
  #       stop(paste("Error during data import: ", filePath, " in line ", as.character(rowNumber), " 
  #                  of the configTable does not correspond to not an existing file!", sep=""))
  #     }
  #     }
}
