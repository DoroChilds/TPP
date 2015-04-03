removeDuplicates = function(inDF, refColName, nonNAColNames, qualColName=NULL){
  ## Remove duplicate entries in ID column
  message("Removing duplicate identifiers using quality column '", qualColName, "'...")
  nonUniques = unique( inDF[duplicated(inDF[refColName]), refColName] )
  
  retDF = subset(inDF, !(get(refColName) %in% nonUniques))
  
  for(nU in nonUniques){
    tmpDF = subset(inDF, get(refColName) == nU)
    nonNArows = NULL
    for(r in 1:nrow(tmpDF)){
      if(any(!is.na(tmpDF[r, nonNAColNames]))){
        nonNArows = c(nonNArows, r)
      }
    }
    if(length(nonNArows) > 1){
      if(is.null(qualColName)){
        useRow = 1
      } else {
        qualVals = tmpDF[nonNArows, qualColName]
        useRow = match(max(qualVals), qualVals)
      }
    } else {
      useRow = nonNArows[1]
    }
    retDF = rbind(retDF, tmpDF[useRow, ])
  }
  message(nrow(retDF), " out of ", nrow(inDF), " rows kept for further analysis.")
  return(retDF)
}
