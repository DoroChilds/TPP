mergeOutputTables_TR <- function(dataList, pValDF, qualCheckDF){
  df1 <- dataList$fcDF # Fold change columns
  df2 <- dataList$curveParDF # Melting curve parameter columns
  df3 <- dataList$plotCol # Plot column
  df4 <- dataList$presenceDF # Which proteins where identified in which experiment?
  df5 <- dataList$modelInfoDF # Additional quality information about the model fits
  df6 <- dataList$otherAnnotDF # Further annotation columns
  
  if (is.null(pValDF)){
    pValDF <- data.frame(Protein_ID = df1$Protein_ID)
  }
  if (is.null(qualCheckDF)){
    qualCheckDF <- data.frame(Protein_ID = df1$Protein_ID)
  }
  
  # Merge all data frames:
  newList <- list(df1, df2, df3, pValDF, qualCheckDF, df4, df5, df6)
  outTable <- join_all(newList, by = "Protein_ID") %>% arrange(Protein_ID)
  message("Results table created successfully.\n")
  return(outTable)
}

mergeOutputTables_CCR <- function(dataList, qualCheckDF){
  ## Initialize output table and append global columns
  ## Initialize data frames to store annotation columns topic wise
  df0 <- dataList$fcOrig
  df1 <- dataList$fcNorm
  df2 <- dataList$fcTransf
  df3 <- dataList$modelPars
  df4 <- dataList$plotCol
  df5 <- dataList$transfDF
  df6 <- dataList$modelInfo
  df7 <- dataList$presenceDF
  df8 <- dataList$otherAnnotDF
  df9 <- dataList$fcRefNorm
  
  ## Merge unmodified/ normalized/ transformed fold-change columns:
  outTable <- df0
  if (!all(is.na(df9[,!grepl("Protein_ID", colnames(df9))]))){
    outTable <- join(outTable, df9, by="Protein_ID")    
  }
  if (!all(is.na(df1[,!grepl("Protein_ID", colnames(df1))]))){
    outTable <- join(outTable, df1, by="Protein_ID")    
  }
  if (!all(is.na(df2[,!grepl("Protein_ID", colnames(df2))]))){
    outTable <- join(outTable, df2, by="Protein_ID")    
  }
  
  ## Add model parameters:
  outTable <- join(outTable, df3, by="Protein_ID")
  
  ## Add plot column:
  outTable <- join(outTable, df4, by="Protein_ID")
  
  ## Add quality check results:
  outTable <- join(outTable, df5, by="Protein_ID")
  
  ## Add information about model fit:
  outTable <- join(outTable, df6, by="Protein_ID")
    
  ## Add columns indicating which proteins where identified in which experiment:
  outTable <- join(outTable, df7, by="Protein_ID")
  
  ## Add futher annotation columns that were directly imported from the input files:
  outTable <- join(outTable, df8, by="Protein_ID")
  
  ## Sort result table alphabetically according to protein ids:
  outTable <- arrange(outTable, Protein_ID)
  message("Results table created successfully.\n")
  
  return(outTable)
}

resultTabFromList <- function(dataList){
  ## to do: migrate construction of the tpptr result table here
  #  nData <- length(dataList)
  outTable <- arrange(join_all(dataList, by="Protein_ID", type="full"), Protein_ID)
  #   if (nData>1){
  #     outTable <- join(dataList[[1]], dataList[[2]], by="Protein_ID")
  #     if (nData>2){
  #       for (n in names(dataList)[3:]){
  #         outTable <- join()
  #       }
  #     }
  #   }
  return(outTable)
}