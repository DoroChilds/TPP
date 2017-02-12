checkAndReturnDataSetting <- function(dataSettings, fieldName, columnNames, sysCall){
  
  ## Check if 'fieldName' exists in the list 'dataSettings' and 
  ## check its content. Particularly, check if it points to an entry in 
  ## 'columnNames'
  
  
  if (!any(fieldName %in% names(dataSettings))){
    
    stop("attr(data, 'importSettings') must contain a field named '", 
         fieldName, "'.")
    
  } else {
    
    fieldEntry <- dataSettings[[fieldName]]
    
  }
  
  charFields <- c("proteinIdCol", "uniqueIdCol", "intensityStr", "nonZeroCols", 
                  "fcStr", "fcStrNorm")
  numFields <- c("r2Cutoff", "fcCutoff", "slopeBounds", "fcTolerance")
  
  if (any(fieldName %in% charFields)){
    
    if (!is.character(fieldEntry)){
      
      stop("attr(data, 'importSettings')$", fieldName, 
           " must be of class character.")
      
    }
  }
  
  if (any(fieldName %in% numFields)){
    
    
    if(!is.numeric(fieldEntry)){
      
      stop("attr(data, 'importSettings')$", fieldName, 
           " must be of class numeric")
      
    }
  }
  
  msgTxt <- paste0(
    "Found the following column name in attr(data, 'importSettings')$", 
    fieldName, ": '", fieldEntry, "'")
  message(msgTxt)
  
  exactMatches <- c("proteinIdCol", "uniqueIdCol", "addCol", "qualColName", 
                    "nonZeroCols")
  substrings <- c("intensityStr", "fcStr", "fcStrNorm")
  
  if (fieldName %in% exactMatches){
    
    if (!all(fieldEntry %in% columnNames)){
      
      stop("The column specified by attr(data, 'importSettings')$", fieldName, 
           " was not found in the column names of 'data'.")
      
    }
    
    
  } else if (fieldName %in% substrings){
    
    if (!any(grepl(fieldEntry, columnNames))){
      stop("The given prefix specified by attr(data, 'importSettings')$", 
           fieldName, " was not found in the column names of 'data'.")
    }
    
  }
  
  return(fieldEntry)
}