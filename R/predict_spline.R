predict_spline <- function(splineModel, x, factorInfo){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  i <- NULL
  
  # Check for missing function arguments if error is not produced by default
  checkFunctionArgs(match.call(), c("splineModel"))
  
  if (length(x) == 0){
    stop("'x' must have at least one element to enable predictions by the smoothing spline model")
  }
  
  if (!is.numeric(x)){
    stop("'x' must be a numeric vector to enable predictions by the smoothing spline model")
  }
  
  
  ## Make sure that newdata contains comparison Factor when
  ## predicting by the alternative model:
  if (inherits(splineModel, "lm")){
    
    if (missing(factorInfo)){
      
      factorInfo <- extract_fit_factors(splineModel = splineModel, mode = "values") 
      
    }
    
    
    if (nrow(factorInfo) > 0){
      
      newDat <- data.frame(x = x, i = 1:length(x)) %>% group_by(i, x) %>% 
        do(factorInfo) %>% ungroup %>% select(-i)
      
    } else {
      
      newDat <- data.frame(x = x)
      
    }
    
    
    ## Start prediction:
    y <- predict(splineModel, newdata = newDat)
    
    out <- newDat %>% mutate(y = y)
    
  } else {
    
    out <- data.frame(x = x, y = NA)
  }
  
  return(out) 
}