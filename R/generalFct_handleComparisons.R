createComparisonTable <- function(infoTable){
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

assignCompNumber_to_expName <- function(compDF, expNames){
  ## Assign a number to each experiment name that indicates to which comparison
  ## pair it belongs. This number will be used to assign matching colors to the
  ## in the experiments in the melting curve plots. Only works when each 
  ## experiment occcurs exactly once in the comparisons. 
  compNumbers <- rep(NA, length(expNames))
  if (!is.null(compDF)){    
    compNames1 <- compDF$testGroup
    compNames2 <- compDF$refGroup
    ## 1. check whether expNames occur in one comparison each
    flag_namesOccurOnlyOnce <- max(table(c(compNames1, compNames2)))==1
    if (flag_namesOccurOnlyOnce){
      ## 2. assign compNumbers according to the comparison a name is involved in
      for (i in 1:nrow(compDF)){
        compNumbers[expNames %in% c(compNames1[i], compNames2[i])] <- i
      }
    }
  }
  return(compNumbers)
}