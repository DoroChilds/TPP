importTR_tidy <- function(configTable, data, idVar, fcStr, naStrs, qualColName, 
                          type){
  message("Importing data...\n")
  
  ## Check configTable for consistency and extract all relevant information:
  configTableContents <- importCheckConfigTable(infoTable=configTable, type=type)
  expInfos <- configTableListToLong(configTableContents)
  expNames   <- configTableContents$expNames
  files      <- configTableContents$files
  compStrs   <- configTableContents$compStrs
  labels     <- configTableContents$labels

  ## If data is given as a list of dataframes, check whether the names are
  ## consistent with the 'Experiment' column in configTable:
  importCheckExperimentNames(expNames = expNames, dataframes = data)
  
  ## Check whether dataframes or filenames are specified for data import:
  argList <- importFct_CheckDataFormat(dataframes = data, files = files, 
                                       expNames = expNames)
  data  <- argList[["dataframes"]]
  files <- argList[["files"]]
  
  if (!is.null(files)){
    ## Import a experiments from files
    data <- importFct_readFiles(files = files, naStrs = naStrs)
  }
  ## Add experiment names and combine to one wide table:
  data2 <- sapply(expNames, simplify = FALSE, USE.NAMES = FALSE,
                  function(nameTmp){
                    return(data[[nameTmp]] %>% mutate(experiment = nameTmp))
                  }) %>%
    bind_rows %>% mutate_if(is.character, factor)
  
  ## Detect fold change column names:
  fcColNames <- importFct_fcCols(datDF = data2, fcPrefix = fcStr, 
                                 labelSuffix = labels)
  
  ## Prepare data frames:
  data3 <- data2 %>% group_by(experiment) %>%
    do({
      nameTmp <- unique(.$experiment)
      message("\nPreparing ",type, " dataset: ", nameTmp)
      fcTidy <- prepareForImport(dataframe = .,
                                 fcColNames = fcColNames,
                                 qualColName  = qualColName,
                                 idVar        = idVar, 
                                 expName      = nameTmp)
    }) %>% mutate(uniqueID = as.factor(uniqueID))
  
  
  ## Extract fold changes:
  valuesPerLabel <- data3 %>% 
    do(importTidyFoldChanges(., fcColNames, type, fcStr)) %>%
    mutate_if(is.character, factor) %>%
    ungroup %>%
    ## Annotate with temperatures/ concentrations per label
    left_join(expInfos, by = c("experiment", "label")) %>%
    arrange(uniqueID, experiment, x)
  # Unit test: are all labelValues (i.e. temperatures) sorted in the output?
  
  ## Extract additional annotation per protein:
  valuesPerProtein <- data3 %>%
    do(importAdditionalCols(., type, fcColNames)) %>%
    ungroup %>%
    arrange(uniqueID, experiment)
  
  
  # ## Store user-specified contrasts in the annotation fields of the newly 
  # ## created ExpressionSets. They will be retrieved later for the statistical
  # ## comparisons:
  # for (n in names(fcListAll)){
  #   annotation(fcListAll[[n]]) <- c(annotation(fcListAll[[n]]), compStrs)     
  # }
  # message("\n")
  
  # hier weiter
  out <- list(concentrations = valuesPerLabel, 
              annotation = valuesPerProtein,
              contrasts = contrasts)
  
  return(fcListAll)
}

configTableListToLong <- function(cfgEntries){
  expNames <- as.character(cfgEntries$expNames)
  labels   <- as.character(cfgEntries$labels)
  conditions  <- cfgEntries$expCond
  conditions2 <- ifelse(is.null(conditions), NA_character_, conditions)
  
  expInfo <- data.frame(experiment = expNames,
                        condition = conditions2) %>%
    mutate_all(as.character)
  
  labelValues <- sapply(1:length(expNames), simplify = FALSE, function(i){
    xTmp <- cfgEntries$tempMatrix[i,]
    data.frame(x = xTmp, label = labels, experiment = expNames[i], 
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  labelInfo <- expand.grid(experiment = cfgEntries$expNames, 
                           label = labels) %>%
    mutate_all(as.character) %>%
    left_join(expInfo, by = "experiment") %>%
    left_join(labelValues, by = c("experiment", "label")) %>%
    mutate_if(is.character, factor) %>%
    arrange(experiment)
  
  return(labelInfo)
}

prepareForImport <- function(dataframe, fcColNames, qualColName, idVar, expName){
  ## Replace NAs in ID variable:
  expName <- as.character(expName)
  dataframe <- data.frame(dataframe) # not a tibble (downstream function still expects scalar vectors when using '$' operator)
  df2 <- importFct_makeUniqueIDs(inDF = dataframe, 
                                 idColName = idVar, 
                                 expName = expName) %>% 
    rename_("uniqueID" = idVar)
  
  ## Make sure that ID variable 'idvar' is unique:
  df3 <- importFct_removeDuplicates(inDF = df2, 
                                    refColName = "uniqueID", 
                                    nonNAColNames = fcColNames, 
                                    qualColName = qualColName)
  
  return(df3)
}

importTidyFoldChanges <- function(dat, fcColNames, type, fcStr){
  nameTmp <- unique(dat$experiment)
  message("Importing fold changes from ",type, " dataset: ", nameTmp)
  fcTidy <- importFct_extractFCs(datDF = dat, colsFC = fcColNames, 
                                 type = NULL) %>%
    as.data.frame  %>%  
    mutate(uniqueID = extract2(dat, "uniqueID")) %>%
    gather_("key", "y", fcColNames) %>%
    mutate(label = gsub(fcStr, "", key),
           columnPrefix = fcStr) 
  importFct_reportValidValues(fcTidy, nameTmp)
  return(fcTidy)
}

importAdditionalCols <- function(dat, type, fcColNames){
  nameTmp <- unique(dat$experiment)
  message("Importing additional row-wise annotation from ",type, " dataset: ", nameTmp)
 ## Retrieve protein annotation columns:
  allCols <- colnames(dat)
  colsAnnot <- setdiff(allCols, fcColNames)
  datAnnot <- dat[, colsAnnot]
  return(datAnnot)
}

importFct_reportValidValues <- function(datLong, expName){
  ## Report size of imported dataset and success rate (how many proteins provide 
  ## enough data to be acutally used for curve fitting?):
  proteinSummary <- datLong %>% group_by(uniqueID) %>% 
    mutate(yIsNotNA = !is.na(y)) %>%
    summarise(sumOfNotNA = sum(yIsNotNA)) %>%
    mutate(atLeast3nonNAs = sumOfNotNA > 2)
  numTotal <- nrow(proteinSummary)
  numValid <- sum(proteinSummary$atLeast3nonNAs)
  message("  -> ", expName, " contains ", numTotal, " proteins.")
  message("  -> ", numValid, " out of ", numTotal, " proteins (",
          round(numValid/numTotal * 100,2), 
          "%) suitable for curve fit (criterion: > 2 valid fold changes per protein).")
  return(NULL)
}