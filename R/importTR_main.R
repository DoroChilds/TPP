importTR_main <- function(configTable, data, idVar, fcStr, naStrs, qualColName, 
                          type){
  ## Wrapper function to invoke the individual steps necessary for data import.
  
  message("Importing data...\n")
  
  ## Check configTable for consistency and extract all relevant information:
  configTableContents <- importCheckConfigTable(infoTable = configTable, 
                                                type = type)
  expNames   <- configTableContents$expNames
  expCond    <- configTableContents$expCond
  files      <- configTableContents$files
  compStrs   <- configTableContents$compStrs
  labels     <- configTableContents$labels
  tempMatrix <- configTableContents$tempMatrix
  
  ## If data is given as a list of dataframes, check whether the names are
  ## consistent with the 'Experiment' column in configTable:
  importCheckExperimentNames(expNames=expNames, dataframes=data)
  
  ## Check whether dataframes or filenames are specified for data import:
  argList <- importFct_CheckDataFormat(dataframes=data, files=files, 
                                       expNames=expNames)
  data  <- argList[["dataframes"]]
  files <- argList[["files"]]
  
  if (!is.null(files)){
    ## Import a experiments from files
    data <- importFct_readFiles(files=files, naStrs=naStrs)
  }
  
  ## Experiment wise pre-processing prior to combination into one object:
  dataFinal <- importFct_preprocessData(data = data, idVar = idVar)
  
  ## Import tables, convert into ExpressionSet format, and store in list:
  fcListAll <- sapply(
    1:length(expNames), simplify=FALSE, USE.NAMES = FALSE,
    function(i){
      importFct_df_to_eSet(dataframe    = dataFinal[[expNames[i]]],
                           labels       = labels,
                           labelValues  = tempMatrix[i,],
                           name         = expNames[i],
                           condition    = expCond[i],
                           idVar        = idVar,
                           fcStr        = fcStr,
                           qualColName  = qualColName,
                           naStrs       = naStrs,
                           type         = type)
    })
  names(fcListAll) <- expNames
  
  ## Store user-specified contrasts in the annotation fields of the newly 
  ## created ExpressionSets. They will be retrieved later for the statistical
  ## comparisons:
  for (n in names(fcListAll)){
    annotation(fcListAll[[n]]) <- c(annotation(fcListAll[[n]]), compStrs)
  }
  message("\n")
  
  return(fcListAll)
}
