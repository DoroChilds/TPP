#' @title Tidy up expressionSets
#' @description Convert list of expressionSets (intermediate output of several 
#' TPP-TR functions) to tidy tables.
#' @param tppESetList A list of expressionSets, returned 
#' by most TPP-TR functions.
#' @param returnType A string with two possible values: "exprs", "featureData".
#' @details expressionSet lists are for example produced by 
#' \code{\link{tpptrImport}}, \code{\link{tpptrNormalize}}, 
#' \code{\link{tpptrCurveFit}}.
#' @return Either the fold changes per protein across all experiments 
#' (if \code{returnType = "exprs"}), or the 
#' additional annotation per protein and experiment (if \code{returnType = "featureData"}). For example, the
#' peptide counts per identified protein can be found here.
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' concentrations <- tpptrTidyUpESets(tpptrData)
#' additionalInfos <- tpptrTidyUpESets(tpptrData, returnType = "featureData")
#' summary(concentrations)
#' @export
tpptrTidyUpESets <- function(tppESetList, returnType = "exprs"){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  experiment = testGroup = refGroup = comparisonFactor = labelValue = 
    foldChange <- NULL
  
  expInfo <- sapply(tppESetList, annotation)
  expNames <- unname(expInfo["name",])
  conditions <- unname(expInfo["condition",])
  
  replicatePattern <- ".*[^[:digit:]]+([:digit:]+)" # Only use the last digits
  
  if (returnType == "exprs"){
    # 1. Create long table of fold changes per protein and TMT-label:
    # 1.1 Convert list of ExpressionSets to long table:
    df1 <- eSetsToLongTable_fc(tppESetList) %>% tibble::as_tibble()
    
    # 1.2 Extract information about condition and replicate from experiment names:
    df2 <- df1 %>% 
      mutate(condition = plyr::mapvalues(experiment, expNames, conditions)) %>%
      extract(experiment, c("replicate"), replicatePattern, remove = FALSE) %>%
      mutate(replicate = paste("Replicate", replicate, sep = "")) # avoid problems due to unnoticed factor-to-numeric conversions
    
    allReplicatesUnique <- df2 %>% distinct(experiment, replicate) %>% 
      extract2("replicate") %>% table() %>% equals(1) %>% all()
    
    if(allReplicatesUnique) {
      df2$replicate = df2$experiment
    }
    
    # ## 1.3 extract info about comparisons to be performed:
    # annot <- sapply(tppESetList, annotation)
    # compRows <- grepl("comparison", rownames(annot), ignore.case=TRUE)
    # if (any(compRows)){
    # compDF <- sapply(tppESetList, annotation) %>%
    #   createComparisonTable() %>% 
    #   select(testGroup, refGroup) %>%
    #   gather(comparisonFactor, experiment, c(testGroup, refGroup)) %>%
    #   mutate(experiment = factor(experiment), 
    #          comparisonFactor = factor(comparisonFactor, 
    #                               levels = c("testGroup", "refGroup")))
    # ## Add column defining the contrasts for the comparisons:
    # df3 <- df2 %>% left_join(compDF, by = "experiment")
    # } else df3 <- df2
    df3 <- df2
    
    # 1.4 Rename some columns to make them more intuitive to understand:
    df4 <- df3  %>% rename(uniqueID = id, x = labelValue, y = foldChange) 
    
    # 1.5 Convert to increase efficiency:
    out <- data.table(df4, key = "uniqueID")
    
  } else if (returnType == "featureData"){
    
    # 2. Create long table with further annotation of each protein:
    # 2.1 Convert list of ExpressionSets to long table:
    df1 <- eSetsToLongTable_fData(tppESetList) %>% tibble::as_tibble() 
    
    # 2.2 Extract information about condition and replicate from experiment names:
    df2 <- df1 %>% 
      mutate(condition = plyr::mapvalues(experiment, expNames, conditions)) %>%
      extract(experiment, c("replicate"), replicatePattern, remove = FALSE) %>%
      mutate(replicate = paste("Replicate", replicate, sep = "")) # avoid problems due to unnoticed factor-to-numeric conversions
    
    # 2.3 Rename some columns to make them more intuitive to understand:
    df3 <- df2  %>% rename(uniqueID = id) 
    
    # 2.4 Convert to increase efficiency:
    out <- data.table(df3, key = "uniqueID") # as.tbl %>% mutate_if(is.character, factor) 
  }
  
  return(out)
}

# out <- list(proteinMeasurements = longTable_measurements,
#             proteinAnnotation = longTable_annotation)
