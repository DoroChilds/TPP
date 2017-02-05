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