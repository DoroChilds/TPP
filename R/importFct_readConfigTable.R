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
        stop("Error during data import: ", cfg, " does not belong to a 
             valid configuration file.")
      }
    } else {
      stop("Error during data import: ", cfg, " does not belong to a 
           valid configuration file.")
    }
    cfg <- tab
  }
  return(cfg)
}

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
}

importFct_CheckDataFormat <- function(files, dataframes, expNames){
  ## Check whether dataframes of filenames are specified for data import
  isDF <- !is.null(dataframes)
  isF  <- !is.null(files)
  isBoth <- isDF & isF
  isNone <- !(isDF | isF)
  if (isBoth){
    stop("Data import function received a filename AND a dataframe object. 
         Please specify only one.")
  } else if (isNone){
    stop("Data import function requires a filename or a dataframe object. 
         Please specify one.")
  }
  
  ## If check was sucessful, continue with data import
  ## Check if import refers to a single table or a list of several tables:
  if (isDF) {
    isClassList <- is.list(dataframes) && !is.data.frame(dataframes)
    isClassDF   <- is.data.frame(dataframes)
    if (isClassList){
      ## If argument "dataframes" is a list, check if it contains only data 
      ## frames. Throw error otherwise.
      classesInList <- unique(sapply(dataframes, class))
      if (length(classesInList)>1 | classesInList[1] !="data.frame"){
        stop("Argument 'dataframes' must be either an object of class 
               'data.frame', or a list of such objects.")
      }
    } else if (isClassDF){
      ## If argument "dataframes" is a single data frame, put it into a list 
      ## for consistency
      dataframes <- list(dataframes)
      names(dataframes) <- expNames
    } else{
      ## If argument "dataframes" is neither a single data frame, nor a list 
      ## of data frames, throw an error.
      stop("Argument 'dataframes' must be either an object of class 
             'data.frame', or a list of such objects.")
    }      
    #     files <- vector("list", length(dataframes))
    #     names(files) <- expNames
  } 
  if (isF)  {
    #     dataframes <- vector("list", length(files))
    #     names(dataframes) <- expNames
    files <- as.character(files)
    names(files) <- expNames
  }
  
  return(list(files=files, dataframes=dataframes))
}

importFct_2Ddataframe <- function(filePath, rowNumber){
  if (file.exists(filePath)){
    ## Determine file extension and import in appropriate format:
    strChunks <- strsplit(filePath, "\\.")[[1]]
    fileExtension <- strChunks[length(strChunks)]
    if (fileExtension=="txt") {
      ## Import table from tab-delimited file
      tab <- read.delim(file=filePath, header=TRUE, check.names=FALSE, as.is=TRUE,
                        stringsAsFactors=FALSE, sep="\t", blank.lines.skip = FALSE)
    } else if (fileExtension=="csv"){
      tab <- read.table(file=filePath, header=TRUE, check.names=FALSE, as.is=TRUE,
                        stringsAsFactors=FALSE, sep=",", blank.lines.skip = FALSE)
    } else if (fileExtension=="xlsx") {
      ## Import table in Excel format
      tab <- read.xlsx(filePath)
    } else {
      stop(paste("Error during data import: ", filePath, " in line ", as.character(rowNumber), " 
                  of the configTable does not correspond to not an existing file!", sep=""))
    }
    return(tab)
  } else{
    stop(paste("Error during data import: ", filePath, " in line ", as.character(rowNumber), " 
                  of the configTable does not correspond to not an existing file!", sep=""))
  }
}