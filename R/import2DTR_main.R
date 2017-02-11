import2DTR_main <- function(configTable, data, idVar, fcStr, addCol, naStrs, 
                            intensityStr, qualColName, nonZeroCols){
  # Creates data frame list for 2D-TPP experiment
  # Returns a list of data frames in which the experimental data are stored 
  # row-wise for each protein.
  
  ## Check whether dataframes or filenames are specified for data import:
  files   <- configTable$Path
  ## Do not consider Path columns with empty strings only.
  if (!is.null(files)){
    if (any(files == "")){
      files <- NULL
    }
  } 
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Experiment = Compound = Temperature = RefCol <- NULL
  
  expNames <- configTable$Experiment
  argList <- importFct_CheckDataFormat(dataframes = data, 
                                       files = files, 
                                       expNames = expNames)
  data  <- argList[["dataframes"]]
  files <- argList[["files"]]
  
  if (!is.null(files)){
    ## Remove redundancies for 2D config tables (2 rows per experiment and filename)
    files2 <- files[!duplicated(names(files))] # files[!duplicated(files)]
    ## Import a experiments from files
    data <- importFct_readFiles(files=files2, naStrs=naStrs)
  }
  
  # loop over all rows of the config table
  configTable %>% group_by(Experiment, Compound, Temperature, RefCol)
  
  iVec <- 1:nrow(configTable)
  dataList <- lapply(iVec, function(iTmp){
    rowTmp <- configTable[iTmp,]
    
    ## Get settings or current experiment
    expTmp <- rowTmp$Experiment
    message("Importing 2D-TPP dataset: ", expTmp)
    
    tTmp <- rowTmp$Temperature     ## Corresponding Temperature value
    
    ## Get corresponding data set
    dataTmp <- data[[expTmp]]
    
    # get columns which correspond to the extracted temperature
    noFCCols <- c("Compound", "Experiment", "Temperature", "RefCol", "Path", "Condition")
    allCols <- colnames(rowTmp)
    labelCols <- setdiff(allCols, noFCCols)
    labelValues <- suppressMessages(rowTmp[,labelCols] %>% as.numeric)
    labelColsNum <- labelCols[!is.na(labelValues)]
    
    # Define all columns that should be imported:
    signalCols <- paste(intensityStr, labelColsNum, sep="")
    relevant.cols <- c(idVar, qualColName, nonZeroCols, addCol, signalCols) %>%
      unique
    if (!is.null(fcStr)){
      fcCols <- paste(fcStr, labelColsNum, sep="")
      relevant.cols <- c(relevant.cols, fcCols)
      dataCols <- fcCols
    } else {
      dataCols <- signalCols
    }
    
    ##  Throw error if any of the relevant column names can not be found in the 
    ## column names of the data frame 
    if (!all(relevant.cols %in% colnames(dataTmp))){
      notFound <- paste(setdiff(relevant.cols, colnames(dataTmp)), collapse = "', '")
      stop("The following columns could not be found: '", notFound, 
           "'. Please check the suffices and the additional column names you have specified.")
    }
    
    ## Experiment wise pre-processing prior to combination into one object:
    dataTmp <- splitIDsIntoSeparateRows(singleDat = dataTmp, idVar = idVar)
    
    ## Make sure that ID variable 'idvar' is unique:
    dataFiltered <- importFct_removeDuplicates(inDF = dataTmp, 
                                               refColName = idVar, 
                                               nonNAColNames = dataCols, 
                                               qualColName = qualColName[1])
    
    # Annotate ids with experiment name
    idsTmp <- as.character(dataFiltered[, idVar])
    idsAnnotated <- paste(expTmp, tTmp, idsTmp, sep="_")
    
    # Select relevant columns and annotate by temperature & experiment
    dataFinal <- dataFiltered %>% subset(select = relevant.cols) %>%
      mutate(temperature = tTmp, experiment = expTmp, unique_ID = idsAnnotated)
    
    return(dataFinal)
  })
  
  ## Create data.list names:
  newNames <- sapply(seq(nrow(configTable)), function(iTmp){
    rowTmp <- configTable[iTmp,]
    # get temperature value
    tTmp <- rowTmp$Temperature
    # get experiment id
    expTmp <- rowTmp$Experiment
    # Define unique experiment name
    newName <- paste(expTmp, tTmp, sep="_")
    return(newName) 
  })
  names(dataList) <- newNames
  
  # remove 0 sumionarea values
  out <- importFct_rmZeroSias(configTable = configTable, 
                              data.list = dataList, 
                              intensityStr = intensityStr)
  
  return(out)  
}
