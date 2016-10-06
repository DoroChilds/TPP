#' @title Tidy up expressionSets
#' @description Convert list of expressionSets (intermediate output of several 
#' TR and CCR functions) to list of tidy tables.
#' @param tppESetList A list of expressionSets, which currently still is the format
#' for intermediate results of most TR and CCR functions.
#' @details expressionSet lists are for example produced by 
#' \code{\link{tpptrImport}}, \code{\link{tpptrNormalize}}, 
#' \code{\link{tpptrCurveFit}}, \code{\link{tppccrImport}}, 
#' \code{\link{tppccrNormalize}}, \code{\link{tppccrCurveFit}}.
#' @return A list with two entries: \code{proteinMeasurements} contains the fold 
#' changes per protein, across all experiments. \code{proteinAnnotation} 
#' contains additional annotation per protein and experiment. For example, the
#' peptide counts per identified protein would go here.
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' longTables <- tpptrData %>% tpptrTidyUpESets
#' concentrations <- longTables %>% extract2("proteinMeasurements")
#' additionalInfos <- longTables %>% extract2("proteinAnnotation")
#' summary(concentrations)
#' @export
tpptrTidyUpESets <- function(tppESetList){
  
  # 1. Create long table of fold changes per protein and TMT-label:
  # 1.1 Convert list of ExpressionSets to long table:
  df1 <- eSetsToLongTable_fc(tppESetList) %>% as.tbl()
  
  # 1.2 Extract information about condition and replicate from experiment names:
  expInfo <- sapply(tppESetList, annotation)
  expNames <- unname(expInfo["name",])
  conditions <- unname(expInfo["condition",])
  replicatePattern <- ".*[^[:digit:]]+([:digit:]+)" # Only use the last digits
  df2 <- df1 %>% 
    mutate(condition = plyr::mapvalues(experiment, expNames, conditions)) %>%
    extract(experiment, c("replicate"), replicatePattern, remove = FALSE) %>%
    mutate(replicate = paste("Replicate", replicate, sep = "")) # avoid problems due to unnoticed factor-to-numeric conversions
  
  ## 1.3 extract info about comparisons to be performed:
  compDF <- sapply(tppESetList, annotation) %>%
    createComparisonTable() %>% 
    select(testGroup, refGroup) %>%
    gather(comparisonFactor, experiment, c(testGroup, refGroup)) %>%
    mutate(experiment = factor(experiment), 
           comparisonFactor = factor(comparisonFactor, 
                                levels = c("testGroup", "refGroup")))
  ## Add column defining the contrasts for the comparisons:
  df3 <- df2 %>% left_join(compDF, by = "experiment")
  
  # 1.4 Rename some columns to make them more intuitive to understand:
  df4 <- df3  %>% rename(uniqueID = id, x = labelValue, y = foldChange) 
  
  # 1.5 Convert character vectors to factors:
  df5 <- as.data.frame(unclass(df4)) %>% as.tbl
  longTable_measurements <-  df5
  
  # 2. Create long table with further annotation of each protein:
  # 2.1 Convert list of ExpressionSets to long table:
  df1 <- eSetsToLongTable_fData(tppESetList) %>% as.tbl() 
  # 2.2 Extract information about condition and replicate from experiment names:
  df2 <- df1 %>% 
    mutate(condition = plyr::mapvalues(experiment, expNames, conditions)) %>%
    extract(experiment, c("replicate"), replicatePattern, remove = FALSE) %>%
    mutate(replicate = paste("Replicate", replicate, sep = "")) # avoid problems due to unnoticed factor-to-numeric conversions
  # 2.3 Rename some columns to make them more intuitive to understand:
  df3 <- df2  %>% rename(uniqueID = id) 
  # 2.4 Convert character vectors to factors:
  df4 <- as.data.frame(unclass(df3)) %>% as.tbl
  longTable_annotation <-  df4
  
  out <- list(proteinMeasurements = longTable_measurements,
              proteinAnnotation = longTable_annotation)
}