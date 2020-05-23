AICc <- function(model){
  
  k <- attr(logLik(model), "df")
  n <- nobs(model)
  aicc <- AIC(model) + 2 * k * (k + 1)/(n - k - 1)
  return(aicc)
  
}
