importCheckConfigTable <- function(infoTable, type){
  ## Check the config table for correctness after import, retrieve different
  ## types of information from it (experiment names, conditions, result path,
  ## comparisons, labels, temperatures), and return them in a list.
  
  ## 1. Check whether obligatory experiment information is provided via data 
  ## frame or spreadsheet file.
  infoTable <- importFct_readConfigTable(cfg=infoTable)
  
  ## 2. Check whether table contains the mandatory column 'Experiment', and 
  ## convert all non-alphanumeric characters to '_' in this column:
  infoTable$Experiment <- importFct_checkExperimentCol(infoTable$Experiment)
  infoTable <- subset(infoTable, Experiment != "") ## todo: test
  
  ## to do: if obsolete 'Replicate' column is given, translate to 'comparison' column
  
  if (type == "TR"){  
    ## 3. If condition column does not exist, assign default values:
    infoTable$Condition <- importFct_checkConditions(infoTable$Condition, nrow(infoTable))    
    
    ## 4. Retrieve user-defined comparisons, check for consistency, and summarize 
    ## them in strings that can be stored in the ExpressionSet annotation fields:
    compStrs <- importFct_checkComparisons(confgTable=infoTable)
    
  } else {
    infoTable$Condition <- NA
    compStrs <- NA
  }
  
  ## 5. Check if table contains mandatory label columns. Stop, if not:
  allCols     <- colnames(infoTable)
  compCols    <- grep("comparison", allCols, value=TRUE, ignore.case=TRUE)
  noLabelCols <- c("Experiment", "Condition", "Path", compCols)
  labels      <- setdiff(allCols, noLabelCols)
  labelStr <- paste(labels, collapse=", ")
  message("The following label columns were detected:\n", labelStr, ".")
  
  ## 6. Retrieve matrix of temperatures to each isobaric label
  temperatures <- infoTable[, labels]
  tempMatrix <- importCheckTemperatures(temp=temperatures)
  
  infoList <- list(expNames = as.character(infoTable$Experiment),
                   expCond  = infoTable$Condition,
                   files    = infoTable$Path,
                   compStrs   = compStrs,
                   labels     = labels,
                   tempMatrix = tempMatrix)
  return(infoList)
}