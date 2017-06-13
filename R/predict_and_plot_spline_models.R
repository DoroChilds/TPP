predict_and_plot_spline_models <- function(dat, fits){
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("dat", "fits")) 
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = testHypothesis <- NULL
  
  ## Predict values across the whole range of the independent variable
  ## (avoids re-fitting by geom_smooth):
  xNew <- seq(min(dat$x), max(dat$x), length.out = 50)
  
  modelPred <- invoke_spline_prediction(fits = fits, x = xNew)
  
  fitFactors <- fits %>% 
    group_by(uniqueID, testHypothesis) %>%
    do({
      out <- tibble()
      if(nrow(.) > 0){
        fitFactors <- extract_fit_factors(splineModel = .$fittedModel[[1]], mode = "names")
        if (length(fitFactors) > 0){
          out <- tibble(factors = fitFactors)
        }}
      out
      }) %>%
    ungroup %>%
    select(-uniqueID) %>%
    distinct
  
  ## Create plot displaying measured and predicted values:
  p <- create_spline_plots(measurements = dat, 
                           predictions = modelPred,
                           colorBy = fitFactors,
                           highlightIDs = c(),
                           highlightTxt = "")
  
  
}
