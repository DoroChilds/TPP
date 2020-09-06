compute_spline_auc <- function(splineModel, xmin, xmax){
  ## Compute area under the spline curve
  
  if (!is.numeric(xmin)) stop("'xmin' must be numeric")
  if (!is.numeric(xmax)) stop("'xmax' must be numeric")
  
  if(length(xmin) != 1) stop("'xmin' must have length 1")
  if(length(xmax) != 1) stop("'xmax' must have length 1")
  
  # Start prediction
  if (!inherits(splineModel, "lm")){
    
    aucTable <- data.frame(auc = NA_real_)
    
  } else {
    
    factorInfo <- extract_fit_factors(splineModel = splineModel, mode = "values")
    factors <- colnames(factorInfo)
    
    if (length(factors) > 0){
      
      if(length(factors) == 1){
        factorInfo <- group_by(factorInfo, !!sym(factors))
      } else {
        factorInfo <- group_by(factorInfo, !!!syms(factors))
      }
      
    } else {
      
      factorInfo <- data.frame()
    }
    
    aucTable <- factorInfo  %>%
      do({
        
        int <- try(stats::integrate(
          
          function(x,m) {
            predict_spline(splineModel = m, x = x, factorInfo = .)$y
          }, 
          
          xmin, xmax, splineModel),
          silent = TRUE)
        
        auc <-  ifelse(test = inherits(int, "try-error"), 
                       yes = NA_real_, 
                       no = int$value)
        
        data.frame(auc = auc)
        
      }) 
  }
  
  return(aucTable)
}

