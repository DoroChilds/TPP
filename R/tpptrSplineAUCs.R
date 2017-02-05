tpptrSplineAUCs <- function(data, fits){
  # Calculate areas under the spline curve (ausc) per model and
  # (if applicable) per factorsH1
  # 
  # Wrapper function for the function 'compute_spline_auc'.
  # 
  # Calls 'compute_spline_auc' sequentially for each row in 'fits' and
  # returns a data frame with all columns from 'fits' except of the
  # fittedModel' column.
  # 
  # The argument 'fits' is the output of 'tpptrFitSplines' and contains the
  # fitted smoothing spline models in the column 'fittedModel'.
  
  # Check for missing function arguments:
  checkFunctionArgs(match.call(), c("data", "fits"))
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = x = successfulFit <- NULL
  
  ## Compute areas under the curve for each row:
  dataGrouped <- data %>% filter(uniqueID %in% fits$uniqueID) %>%
    group_by(uniqueID) 
  
  xRanges <- dataGrouped %>% 
    dplyr::summarise(xMin = min(x, na.rm = TRUE), xMax = max(x, na.rm = TRUE))
  
  fits2 <- fits %>% 
    left_join(xRanges, by = "uniqueID")
  
  if (! "fittedModel" %in% colnames(fits)){
    stop("'fits' must contain a column called 'fittedModel'")
  }
  
  message("Calculating areas under the curves for the null and alternative models.")
  aucTable <- fits2 %>%
    filter(successfulFit) %>%
    rowwise() %>% # each row contains one model and will be processed sequentially
    do({
      rowContents <- .
      ## Predict values:
      res <- compute_spline_auc(splineModel = rowContents$fittedModel, 
                                xmin = rowContents$xMin, xmax = rowContents$xMax)
      # Add information about grouping variables. This information was lost
      # when invoking the prediction after rowwise grouping:
      otherCols <- rowContents %>% inset2("fittedModel", NULL) %>% data.frame()
      out <- bind_cols(res, mefa:::rep.data.frame(otherCols, nrow(res)))   
    }) %>% 
    ungroup
  
  # replace NAs in color column by entry 'null model' for null models
  allModels <- fits2 %>% extract2("fittedModel")
  fitFactors <- allModels %>% 
    purrr::map(function(m) {
      if (inherits(m, "lm")){
        return(extract_fit_factors(m, mode = "names"))
        #return(colnames(extract_fit_factors(m)))
      } else {
        return(c())
      }
    }) %>%
    unlist %>% unique
  
  isNullModel <- aucTable$testHypothesis == "null"
  for (factorTmp in fitFactors){
    isNAValue <- is.na(aucTable[[factorTmp]])
    aucTable[which(isNAValue & isNullModel), factorTmp] <- "null model"
  }
  
  return(aucTable)
}

