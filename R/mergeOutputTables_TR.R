mergeOutputTables_TR <- function(dataList, pValDF, qualCheckDF){
  ## Generate final TR output table.
  ## Concatenate individual data frames produced by upstream functions
  ## into a wide table that can be exported to Excel.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Protein_ID <- NULL
  
  ## Check for missing function arguments
  checkFunctionArgs(match.call(), c("dataList"))
  
  ## Check if all relevant fields exist in 'dataList':
  expectedCols <- c("fcDF", "curveParDF", "plotCol", "presenceDF", "modelInfoDF", "otherAnnotDF")
  
  colExists <- expectedCols %in% names(dataList)
  
  if (!all(colExists)){
    
    stop("The following fileds are missing in dataList: '", 
         paste(expectedCols[!colExists], collapse = "', '"), "'")
    
  }
  
  allIDs <- dataList$fcDF$Protein_ID
  
  if (is.null(pValDF)){
    pValDF <- data.frame(Protein_ID = allIDs)
  }
  
  if (is.null(qualCheckDF)){
    qualCheckDF <- data.frame(Protein_ID = allIDs)
  }
  
  ## Retrieve individual data frames and merge them in the desired order
  df1 <- dataList$fcDF # Fold change columns
  df2 <- dataList$curveParDF # Melting curve parameter columns
  df3 <- dataList$plotCol # Plot column
  df4 <- dataList$presenceDF # Which proteins where identified in which experiment?
  df5 <- dataList$modelInfoDF # Additional quality information about the model fits
  df6 <- dataList$otherAnnotDF # Further annotation columns (appended last)
  
  
  # Merge all data frames:
  newList <- list(df1, df2, df3, pValDF, qualCheckDF, df4, df5, df6)
  
  outTable <- join_all(newList, by = "Protein_ID") %>% arrange(Protein_ID)
  
  message("Results table created successfully.\n")
  
  return(outTable)
}