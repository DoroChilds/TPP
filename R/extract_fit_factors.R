extract_fit_factors <- function(splineModel, mode){
  
  ## Define empty return values. Will be returned if not overwritten afterwards.
  ## Check mode argument for correctnes at the same time.
  if (mode == "names"){
    out <- character()
  } else if (mode == "values"){
    out <- data.frame()
  } else {
      stop("argument mode must be one of the following strings: 'names', 'values'")
  }
  
  ## Check whether spline model is a linear model fit:
  if (inherits(splineModel, "lm")){
    
    factorColumns <- splineModel$model %>% 
      dplyr::select(contains("factor"))
    
    colnames(factorColumns) <- colnames(factorColumns) %>% 
      gsub("factor\\(", "", .) %>%
      gsub("\\)", "", .)
    
    
    if (mode == "names"){
      
      out <- colnames(factorColumns)
      
    } else if (mode == "values"){
      
      if (ncol(factorColumns) > 0){
        
        out <- distinct(factorColumns) %>%
          mutate_all(as.character)
        
      }
    } 
  }
  
  return(out)
}
