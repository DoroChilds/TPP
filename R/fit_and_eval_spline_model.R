fit_and_eval_spline_model <- function(df, xNew, modelFormula){
  splineFit <- try(rlm(as.formula(modelFormula), 
                       data = df, maxit = 150),
                   silent = TRUE)
  if (!inherits(splineFit, "try-error")){
    yNew <- predict(splineFit, list(x = xNew))
  } else {
    yNew <- NA
  }
  normDF <- data.frame(x = xNew, splinePrediction = yNew)
  return(normDF)
}