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
    
    m_augmented <- splineModel %>% broom::augment()
    
    factorColumns <- m_augmented %>% 
      colnames %>% 
      grep("factor.", ., value = TRUE) 
    
    factorNames <- factorColumns %>% 
      gsub("factor", "", .) %>% 
      gsub("\\.", "", .)
    
    if (mode == "names"){
      
      out <- factorNames
      
    } else if (mode == "values"){
      
      if (length(factorNames) > 0){
        
        factorValues <- m_augmented %>% 
          subset(select = factorColumns) %>% 
          distinct %>%
          set_names(factorNames) %>%
          mutate_all(as.character)
        
        out <- factorValues
        
      }
    } 
  }
  
  return(out)
}
