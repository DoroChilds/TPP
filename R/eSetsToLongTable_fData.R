eSetsToLongTable_fData <- function(data){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  experiment = variable <- NULL
  
  expNames <- names(data)
  longTabAnnot <- c()
  for (en in expNames){
    datTmp <- data[[en]]
    
    # Retrieve protein ids, pData, fData and fold changes:
    ids <- featureNames(datTmp)
    fDatWide <- pData(featureData(datTmp))
    
    # Long table of fold changes:
    colnames(fDatWide) <- gsub("([^[:alnum:]])", "_", colnames(fDatWide))
    fDatLong <- fDatWide %>%
      mutate(id = ids) %>%
      mutate_all(as.character) %>%
      gather_("variable", "value", colnames(fDatWide)) %>%
      mutate(experiment = en)
    longTabAnnot <- rbind(longTabAnnot, fDatLong)
  }
  
  longTabAnnot <- longTabAnnot %>% 
    arrange(id) %>%
    mutate(id = factor(id),
           experiment = factor(experiment),
           variable = factor(variable))
  
  return(longTabAnnot)
}
