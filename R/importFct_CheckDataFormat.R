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
      classesInList <- dataframes %>% sapply(. %>% inherits(., "data.frame"))
      if (!all(classesInList)){
        stop("Argument 'dataframes' contains elements that are not of type 'data.frame' at the following positions: ", 
                which(!classesInList) %>% paste(collapse = ", "), ".")
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
    #     stop("Please hand over the experimental data as a list of data frames or specifiy respective 
    # file paths in the config table!")
  }
  
  return(list(files=files, dataframes=dataframes))
  }
