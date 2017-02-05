mergeOutputTables_CCR <- function(dataList, qualCheckDF){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Protein_ID <- NULL
  
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