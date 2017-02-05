splitTidyMeasurementsForExport <- function(measurements, proteinInfos){
  # Convert the data obtained by the 'tpptrSplineFitAnd Test' function into a 
  # format that can be used by the 'mergeOutputTables_TR' function in order
  # to create the final result table.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = experiment = colName = newColName = y = variable = value <- NULL
  
  emptyDF <- data.frame(Protein_ID = measurements$uniqueID) %>% distinct
  
  # Convert measurements from long to wide:
  fcDF <- measurements %>% 
    arrange(uniqueID, experiment, colName) %>% # to do: replace the 'colName' column by a more meaningfull name (created by function 'eSetsToLongTable_fc')
    mutate(newColName = paste(colName, experiment, sep = "_")) %>%
    mutate(newColName = factor(newColName, levels = unique(newColName))) %>%
    select(uniqueID, y, newColName) %>% 
    spread(newColName, y) %>% 
    rename(Protein_ID = uniqueID)
  
  # Convert further protein annotation from long to wide:
  if (!is.null(proteinInfos)){
    rmCols <- c("plot", 
                meltCurveParamNames(returnParNames = TRUE, 
                                    returnPerformanceInfo = TRUE))
    
    infoTabfiltered <- proteinInfos %>% 
      arrange(uniqueID, experiment, variable) %>%
      filter(!variable %in% rmCols)
    
    if (nrow(infoTabfiltered) > 0){
      otherAnnotDF  <- infoTabfiltered %>%
        mutate(newColName = paste(variable, experiment, sep = "_")) %>%
        mutate(newColName = factor(newColName, levels = unique(newColName))) %>%
        select(uniqueID, value, newColName) %>% 
        spread(newColName, value) %>% 
        rename(Protein_ID = uniqueID)
    } else {
      otherAnnotDF <- emptyDF
    }
  } else {
    otherAnnotDF <- emptyDF
  }
  
  ## Store in a list that can be used by the 'mergeOutputTables_TR' function.
  dataList <- list(fcDF = fcDF,
                   curveParDF = emptyDF,
                   plotCol = emptyDF,
                   presenceDF = emptyDF,
                   modelInfoDF = emptyDF,
                   otherAnnotDF = otherAnnotDF)
  
  return(dataList)
}
