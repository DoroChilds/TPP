createComparisonTable <- function(infoTable){
  ## Create a table that summarizes the experiment-wise comparisons
  ## specified by the user.
  
  compRows <- grepl("comparison", rownames(infoTable), ignore.case=TRUE)
  
  if (any(compRows)){
    compStrs <- unique(c(infoTable[compRows,]))
    compGroups <- as.data.frame(strsplit(compStrs, "_vs_"), 
                                row.names=c("testGroup", "refGroup"), 
                                stringsAsFactors = FALSE)
    compGroups <- as.data.frame(t(compGroups), row.names=1:ncol(compGroups), 
                                stringsAsFactors = FALSE)
    compGroups$name <- compStrs
    
    ## Check whether the experiments in a pair were also defined as a Vehicle and 
    ## Treatment groups by the user. This information is important for the quality 
    ## checks before and after p-value calculation:
    refConditions  <- infoTable["condition",compGroups$refGroup]
    testConditions <- infoTable["condition",compGroups$testGroup]
    compGroups$refIsVehicle    <- ifelse(is.na(refConditions), FALSE, refConditions == "Vehicle")
    compGroups$testIsTreatment <- ifelse(is.na(testConditions), FALSE, testConditions == "Treatment")
    compGroups$testIsVehicle   <- ifelse(is.na(testConditions), FALSE, testConditions == "Vehicle")
  } else{
    compGroups <- NULL
  }
  return(compGroups)
}