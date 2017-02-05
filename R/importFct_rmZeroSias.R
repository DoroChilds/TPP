importFct_rmZeroSias <- function(configTable, data.list, intensityStr){
  ## @title Remove rows with zero sumionarea values 
  ## 
  ## @description Removes zero sumionare values in a specified data.list so that no errors are
  ##   generated in the following fold change computation step. A corresponding data.list  
  ##   with NAs instead of zeros is returned.
  ##   
  ## @return A list of data frames with NAs instead of zeros.
  ##
  ## @param configTable data frame that specifies important details of the 2D-TPP experiment.
  ## @param data.list list of data frames of corresponding experiment data  
  ## @param intensityStr character string indicating which columns contain the sumionarea 
  ##   values. Those column names containing the suffix \code{intensityStr} 
  ##   will be regarded as containing sumionarea values.  
  out <- lapply(names(data.list), function(l.name){
    
    datTmp <- data.list[[l.name]]
    colsTmp <- colnames(datTmp)
    
    intensity.cols <- grep(intensityStr, colsTmp, value = TRUE)
    
    intensity.df <- subset(datTmp, select = intensity.cols) %>%
      mutate_all(as.character) %>% mutate_all(as.numeric) 
    
    new.intensity.df <- intensity.df %>% mutate_all(replaceZeros)
    
    datTmp[,intensity.cols] <- new.intensity.df 
    
    return(datTmp)
  })
  names(out) <- names(data.list)
  return(out)
}

replaceZeros <- function(x){
  x[which(x == 0)] <- NA
  return(x)
}
