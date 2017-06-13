modelSelector <- function(fitStats, criterion, hypothesis) {
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = splineDF <- NULL
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("fitStats", "criterion", "hypothesis"))
  
  # Check if all relevant colnames exist:
  expectedCols <- c("uniqueID", criterion, "testHypothesis", "splineDF")
  
  colExists <- expectedCols %in% colnames(fitStats)
  
  if (!all(colExists)){
    
    stop("The following columns are missing in fitStats: '", 
         paste(expectedCols[!colExists], collapse = "', '"), "'")
    
  }
  
  # Check column format and contents:
  if (!is.numeric(fitStats$splineDF)) 
    stop("'fitStats$splineDF' must be numeric")
  
  if (!is.numeric(fitStats[[criterion]])) 
    stop("'fitStats[[",criterion,"]]' must be numeric")
  
  if (any(is.na(fitStats$testHypothesis))) 
    warning("'fitsStats$testHypothesis' contains missing values")
  
  if (! hypothesis %in% fitStats$testHypothesis) 
    stop("Given argument '", hypothesis, 
         "' not found in column 'fitsStats$testHypothesis'")
  
  # Replace NA entries in filter column by numerics so that protein is not discarded :
  fitStats <- ungroup(fitStats)
  allIDs <- distinct(fitStats, uniqueID)
  
  out <- fitStats %>% 
    filter_(paste0('testHypothesis == "', hypothesis, '"')) %>%
    mutate(splineDF = ifelse(is.na(splineDF), Inf, splineDF)) %>%
    group_by(uniqueID) %>% 
    filter_(paste0(criterion, ' == min(', criterion, ', na.rm = TRUE)')) %>% 
    dplyr::summarize(splineDF = min(splineDF)) %>% # in case of ties, use the least complex model
    mutate(splineDF = ifelse(is.infinite(splineDF), NA_real_, splineDF)) %>%
    arrange(uniqueID) %>%
    ungroup() %>%
    select(uniqueID, splineDF) %>%
    right_join(allIDs, by = "uniqueID") # Join back proteins with NA in all criterion values. Will receive splineDF = NA
  
  return(out)
}


