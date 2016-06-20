mergeOutputTables_TR <- function(dataList, pValDF, qualCheckDF){
  ## Generate final TR output table.
  
  df1 <- dataList$fcDF
  df2 <- dataList$curveParDF
  df3 <- dataList$plotCol
  df4 <- dataList$presenceDF
  df5 <- dataList$modelInfoDF
  df6 <- dataList$otherAnnotDF
  
  ## Merge fold change column and melting curve parameter columns:
  outTable <- join(df1, df2, by="Protein_ID")
  
  ## Add plot column:
  outTable <- join(outTable, df3, by="Protein_ID")
  
  ## Add columns for p-values and quality checks:
  if (!is.null(pValDF)){
    outTable <- join(outTable, pValDF, by="Protein_ID")
  }
  if (!is.null(qualCheckDF)){
    outTable <- join(outTable, qualCheckDF, by="Protein_ID")
  }
  
  ## Add columns indicating which proteins where identified in which experiment:
  outTable <- join(outTable, df4, by="Protein_ID")
  
  ## Add columns with additional quality information about model fit (boolean
  ## variables indicating whether sufficient non-missing values were
  ## available for model fit and whether the model converged successfully)
  outTable <- join(outTable, df5, by="Protein_ID")
    
  ## Add futher annotation columns that were directly imported from the input files:
  outTable <- join(outTable, df6, by="Protein_ID")
  
#   ## Replace NAs by FALSE in boolean columns:
#   outTable <- resultTab_fct_removeNAsFromBooleanCols(outTable)
  
  ## Sort result table alphabetically according to protein ids:
  outTable <- arrange(outTable, outTable$Protein_ID)
  message("Results table created successfully.\n")
  
  return(outTable)
}

mergeOutputTables_CCR <- function(dataList, qualCheckDF){
  ## Generate final CCR output table.
  
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
  
#   ## Replace NAs by FALSE in boolean columns:
#   outTable <- resultTab_fct_removeNAsFromBooleanCols(outTable)
  
  ## Sort result table alphabetically according to protein ids:
  outTable <- arrange(outTable, outTable$Protein_ID)
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

# resultTab_fct_removeNAsFromBooleanCols <- function(tab){
#   ## Replace NAs by FALSE in boolean columns:
#   
#   ## detect boolean columns:
#   boolPos   <- sapply(tab, is.logical, simplify = TRUE)
#   boolNames <- colnames(tab)[boolPos]
#   for (bn in boolNames){
#     xOld = xNew <- tab[,bn]
#     naPos <- is.na(xOld)
#     xNew[naPos] <- FALSE
#     tab[,bn] <- xNew
#   }
#   return(tab)
# }