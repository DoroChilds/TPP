invoke_spline_prediction <- function(fits, x){
  # Wrapper function for the function 'predict_spline'.
  # 
  # Calls 'predict_spline' sequentially for each row in 'fits' and
  # returns a data frame with all columns from 'fits' except of the
  # fittedModel' column.
  # 
  # The argument 'fits' is the output of 'tpptrFitSplines' and contains the
  # fitted smoothing spline models in the column 'fittedModel'.
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("fits", "x")) 
  
  if (! "fittedModel" %in% colnames(fits)){
    stop("'fits' must contain a column called 'fittedModel'")
  }
  
  predictions <- fits %>%
    rowwise() %>% # each row contains one model and will be processed sequentially
    do({
      rowContents <- .
      ## Predict values:
      res <- predict_spline(splineModel = rowContents$fittedModel, x = x)
      # Add information about grouping variables. This information was lost
      # when invoking the prediction after rowwise grouping:
      otherCols <- rowContents %>% inset2("fittedModel", NULL) %>% data.frame()
      out <- cbind(res, otherCols)
    }) %>% 
    ungroup
  
  # replace NAs in color column by entry 'nullModel' for null models
  allModels <- fits %>% extract2("fittedModel")
  fitFactors <- allModels %>% 
    purrr::map(function(m) {
      if (inherits(m, "lm")){
        return(extract_fit_factors(m, mode = "names")) #return(colnames(extract_fit_factors(m)))
      } else {
        return(c())
      }
    }) %>%
    unlist %>% unique
  
  isNullModel <- predictions$testHypothesis == "null"
  for (factorTmp in fitFactors){
    isNAValue <- is.na(predictions[[factorTmp]])
    predictions[which(isNAValue & isNullModel), factorTmp] <- "null model"
  }
  
  return(predictions)
}