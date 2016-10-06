importCheckConfigTable <- function(infoTable, type){
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
    
  } else if (type == "2D"){
    ## 4. Retrieve user-defined comparisons, check for consistency, and summarize 
    ## them in strings that can be stored in the ExpressionSet annotation fields:
    compStrs <- NA
    
  } else {
    infoTable$Condition <- NA
    compStrs <- NA
  }
  
  if (type == "2D"){ 
    ## 5. Check if table contains mandatory label columns. Stop, if not:
    allCols     <- colnames(infoTable)
    noLabelCols <- c("Compound", "Experiment", "Temperature", "RefCol", "Path")
    labels      <- setdiff(allCols, noLabelCols)
    labelStr <- paste(labels, collapse=", ")
    message("The following label columns were detected:\n", labelStr, ".")
    temperatures <- infoTable$Temperature
    # stop if temperature indication not sufficient for analysis
    if (is.null(temperatures) | length(temperatures)<2){
      stop("2D-TPP analysis cannot be performed with insufficient temperature data values!\n Check 
           whether your configuration table has the correct column names!")
    } else if (length(which(!infoTable$RefCol %in% labels))!=0){
      stop("Each reference column must among the label columns!")
    } else{
      return(infoTable)
    }
    
  } else {    
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
}