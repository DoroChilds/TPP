importFct_makeUniqueIDs = function(inDF, idColName, expName){
  ## Make sure that a column with the name stored in idVar exists:
  if (!any(colnames(inDF) == idColName)){
    stop("idVar Column '", idColName, "' not found in dataset '", expName, "'.")
  }
  
  idColumn <- as.character(inDF[[idColName]])
  countSuffix <- seq(1:length(idColumn[is.na(idColumn)]))
  idColumn[is.na(idColumn)] <- paste('NA_in', expName, countSuffix, sep='_')
  inDF[[idColName]] <- idColumn
  return(inDF)
}