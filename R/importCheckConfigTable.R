importCheckConfigTable <- function(infoTable, type){
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("infoTable", "type"))
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Experiment = Path = Compound <- NULL
  
  # Check input arguments and produce informative error messages:
  isValidDF <- FALSE
  if(is.data.frame(infoTable)){
    if ((ncol(infoTable) > 1) & ("Experiment" %in% colnames(infoTable))){
      isValidDF <- TRUE
    }
  }
  
  if (!is.character(infoTable) & !isValidDF){
    stop("'infoTable' must either be a data frame with an 'Experiment' column 
         and at least one isobaric label column, or a filename pointing at a 
         table that fulfills the same criteria")
  }
  
  isValidType <- type %in% c("TR", "CCR", "2D")
  if (!isValidType) {
    stop("'type' must have one of these values: 'TR', 'CCR', '2D'")
  }
  
  
  ## 1. Check whether obligatory experiment information is provided via data 
  ## frame or spreadsheet file.
  infoTable <- importFct_readConfigTable(cfg=infoTable)
  
  ## 2. Check whether table contains the mandatory column 'Experiment', and 
  ## convert all non-alphanumeric characters to '_' in this column:
  infoTable$Experiment <- importFct_checkExperimentCol(infoTable$Experiment)
  infoTable <- subset(infoTable, Experiment != "")
  
  ## 2.1 Check whether table contains a column 'Path', and remove if empty:
  givenPaths <- NULL
  if (any("Path" %in% colnames(infoTable))) {
    if (all(infoTable$Path == "") || all(is.na(infoTable$Path))){
      message("Removing empty 'Path' column from config table")
      infoTable <- infoTable %>% select(-Path)
    } else {
      givenPaths <- infoTable$Path
    }
  }
  
  ## 3. Retrieve user-defined comparisons, check for consistency, and summarize 
  ## them in strings that can be stored in the ExpressionSet annotation fields
  ## (TR part only):
  if (type == "TR"){  
    compStrs <- importFct_checkComparisons(confgTable=infoTable)
  } else {
    compStrs <- NA
  }
  
  ## 4. If condition column does not exist, assign default values 
  ## (TR part only):
  if (type == "TR"){  
    infoTable$Condition <- importFct_checkConditions(infoTable$Condition, nrow(infoTable))    
  }  else {
    infoTable$Condition <- NULL
  }
  
  
  ## 5. Check if table contains mandatory label columns. Stop, if not:
  allCols     <- colnames(infoTable)
  
  labelCols <- detectLabelColumnsInConfigTable(allColumns = allCols)
  
  ## 6. Remove label columns that do not contain at least 1 number:
  labelValues <- infoTable[,labelCols]
  labelValuesNum <- suppressWarnings(labelValues %>% apply(2, as.numeric))
  if (is.matrix(labelValuesNum)) {
    isInvalid <- labelValuesNum %>% apply(2, is.na) %>% apply(2, all)
  } else if (is.vector(labelValuesNum)){
    isInvalid <- is.na(labelValuesNum)
  }
  invalidLabels <- labelCols[isInvalid]
  infoTable[,invalidLabels] <- NULL
  
  labelColsNew <- labelCols[!isInvalid]
  # infoTable[,labelColsNew] <- labelValuesNum[,labelColsNew]
  
  labelStr <- paste(labelColsNew, collapse=", ")
  message("The following valid label columns were detected:\n", labelStr, ".")
  
  if (type == "2D"){ 
    ## 2D TPP specific checks:
    ## 7. Check temperature values (TR-part and 2D-TPP part):
    temperatures <- infoTable$Temperature
    # stop if temperature indication not sufficient for analysis
    if (is.null(temperatures) | length(temperatures)<2){
      m1 <- "Insufficient temperatures (<2) specified in config file." 
      m2 <- "Does your configuration table have the correct column names?"
      stop(m1, "\n", m2)
    } else if (length(which(!infoTable$RefCol %in% labelColsNew))!=0){
      stop("Labels in reference column not found in any of teh label columns.")
    } 
    ## 8. Check if table contains mandatory column "Compound" and remove special 
    ##    characters:
    hasCompoundCol <- any(allCols == "Compound")
    if (!hasCompoundCol){
      m <- "Config table of a 2D-TPP experiment needs a 'Compound' column."
      stop(m, "\n")
    } else {
      infoTable <- infoTable %>% 
        mutate(Compound = gsub("([^[:alnum:]])", "_", Compound))
    }
    out <- infoTable
  } else {   
    ## TR- and CCR- specific checks:
    ## 6. Retrieve matrix of temperatures to each isobaric label
    temperatures <- subset(infoTable, select = labelColsNew)
    tempMatrix <- importCheckTemperatures(temp=temperatures)
    
    infoList <- list(expNames = as.character(infoTable$Experiment),
                     expCond  = infoTable$Condition,
                     files    = givenPaths,
                     compStrs   = compStrs,
                     labels     = labelColsNew,
                     tempMatrix = tempMatrix)
    out <- infoList
  }
  return(out)
}
