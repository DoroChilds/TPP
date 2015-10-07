importCheckConfigTable <- function(infoTable, type){
  ## 1. Check whether obligatory experiment information is provided via data 
  ## frame or spreadsheet file.
  infoTable <- importFct_readConfigTable(cfg=infoTable)
  
  ## 2. Check whether table contains the mandatory column 'Experiment':  
  if (!"Experiment" %in% names(infoTable)){
    m<-"Config table needs an 'Experiment' column with unique experiment IDs."
    stop(m, "\n")
  }
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

importCheckExperimentNames <- function(expNames, dataframes){
  ## If data is given as a list of dataframes, check whether the names are
  ## consistent with the 'Experiment' column in configTable (argument expNames):
  if (is.list(dataframes) && !is.data.frame(dataframes)){
    if (!identical(expNames, names(dataframes))){
      stop("The names of the data objects in 'data' differ from the names given in the Experiment column of 'configTable'.")
    }
  }
}

importCheckTemperatures <- function(temp){
  tempMatrix <- as.matrix(temp)
  ## Make sure that temperature matrix has non-missing values:
  naRows <- apply(is.na(tempMatrix), 1, all)
  
  if(any(naRows)){
    stop("Row(s) ", paste(which(naRows), collapse=", "), " in the configuration table contain only missing temperature values.")
  }
  
  return(tempMatrix)
}

importFct_checkComparisons <- function(confgTable){
  ## Check, if comparisons were specified by the user (via the 
  ## 'Comparison' column in the config table). IF yes, check them for 
  ## consistency, and summarize them in strings that can be stored in 
  ## the ExpressionSet annotation fields.
  
  ## Preparation:
  expConds <- confgTable$Condition
  expNames <- confgTable$Experiment
  
  ## 1. Check whether table contains 'Comparison' columns:
  compCols <- grep("Comparison", colnames(confgTable), ignore.case=TRUE, 
                   value=TRUE)
  
  ## 2. Check whether these columns specify valid comparisons (exactly two 
  ## alphanumeric values per column):
  compChars <- apply(confgTable[compCols], 2, function(x){
    length(grep("[[:alnum:]]", x, value=TRUE))})
  
  
  ## 3. Produce a warning, if a column with the prefix 'comparison' does not 
  ## contain exactly two entries:
  comp_unequal_two <- compChars != 2
  if (any(comp_unequal_two)){
    warning("\nThe following comparison columns could not be evaluated because they did not contain exactly two entries:\n\t\t", 
            paste(compCols[comp_unequal_two], collapse=",\n\t\t"))
  }
  
  ## 4. Create characters that describe the comparisons to be performed and that 
  ## can be stored in the ExpressionSet annotation field:
  validCompCols <- compCols[!comp_unequal_two]
  
  allCompStrs <- c()  
  if (length(validCompCols) > 0){
    message("Comparisons will be performed between the following experiments:")
    for (colName in validCompCols){
      current_compEntries <- confgTable[,colName]
      current_compRows    <- grep("[[:alnum:]]", current_compEntries)
      current_compExps    <- expNames[current_compRows]
      compRef    <- current_compExps[1]
      compTreatm <- current_compExps[2]
      
      ## Re-define using condition info (if avaliable):
      if ("Condition" %in% names(confgTable)){
        current_compConds <- expConds[current_compRows]
        if ("Vehicle" %in% current_compConds && "Treatment" %in% current_compConds){
          compRef    <- current_compExps[current_compConds == "Vehicle"]
          compTreatm <- current_compExps[current_compConds == "Treatment"]
        }
      }
      
      ## Store and report comparison:
      compStr <- paste(compTreatm, "_vs_", compRef, sep="")
      names(compStr) <- colName
      message(compStr)
      allCompStrs <- c(allCompStrs, compStr)
    }
    message("\n")
  }
  
  return(allCompStrs)
}

importFct_checkConditions <- function(condInfo, expectedLength){
  ## Assign generic condition information to each experiment, if not specified 
  ## by the user.
  flagGenerateConds <- FALSE
  if (is.null(condInfo)){
    message("No information about experimental conditions given.")
    flagGenerateConds <- TRUE
  } else{
    condInfo <- as.character(condInfo)
    condLevels <- unique(condInfo)
    if (!identical(sort(condLevels), c("Treatment", "Vehicle"))){
      message("Condition column contains the values '", paste(sort(condLevels), collapse="', '"), "'.")
      message("Need the values 'Treatment' and 'Vehicle' instead.")
      flagGenerateConds <- TRUE
    }
  }
  if (flagGenerateConds){
    message("Assigning NA to all conditions.
Reminder: recognition of Vehicle and Treatment groups during pairwise 
comparisons is only possible when they are specified in the config table.\n")
    condInfo <- rep(NA_character_, expectedLength)
  }
  return(as.character(condInfo))
}

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