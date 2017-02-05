importFct_checkComparisons <- function(confgTable){
  ## Check, if comparisons were specified by the user (via the 
  ## 'Comparison' column in the config table). IF yes, check them for 
  ## consistency, and summarize them in strings that can be stored in 
  ## the ExpressionSet annotation fields.
  
  ## Preparation:
  expConds <- confgTable$Condition
  expNames <- confgTable$Experiment
  
  ## 1. Check whether table contains 'Comparison' columns:
  compCols <- grep("Comparison", colnames(confgTable), ignore.case=TRUE, 
                   value=TRUE)
  
  ## 2. Check whether these columns specify valid comparisons (exactly two 
  ## alphanumeric values per column):
  compChars <- apply(confgTable[compCols], 2, function(x){
    length(grep("[[:alnum:]]", x, value=TRUE))})
  
  
  ## 3. Produce a warning, if a column with the prefix 'comparison' does not 
  ## contain exactly two entries:
  comp_unequal_two <- compChars != 2
  if (any(comp_unequal_two)){
    warning("\nThe following comparison columns could not be evaluated because they did not contain exactly two entries:\n\t\t", 
            paste(compCols[comp_unequal_two], collapse=",\n\t\t"))
  }
  
  ## 4. Create characters that describe the comparisons to be performed and that 
  ## can be stored in the ExpressionSet annotation field:
  validCompCols <- compCols[!comp_unequal_two]
  
  allCompStrs <- c()  
  if (length(validCompCols) > 0){
    message("Comparisons will be performed between the following experiments:")
    for (colName in validCompCols){
      current_compEntries <- confgTable[[colName]]
      current_compRows    <- grep("[[:alnum:]]", current_compEntries)
      current_compExps    <- expNames[current_compRows]
      compRef    <- current_compExps[1]
      compTreatm <- current_compExps[2]
      
      ## Re-define using condition info (if avaliable):
      if ("Condition" %in% names(confgTable)){
        current_compConds <- expConds[current_compRows]
        if ("Vehicle" %in% current_compConds && "Treatment" %in% current_compConds){
          compRef    <- current_compExps[current_compConds == "Vehicle"]
          compTreatm <- current_compExps[current_compConds == "Treatment"]
        }
      }
      
      ## Store and report comparison:
      compStr <- paste(compTreatm, "_vs_", compRef, sep="")
      names(compStr) <- colName
      message(compStr)
      allCompStrs <- c(allCompStrs, compStr)
    }
    message("\n")
  }
  
  return(allCompStrs)
}
