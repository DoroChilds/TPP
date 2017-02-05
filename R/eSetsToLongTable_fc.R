eSetsToLongTable_fc <- function(data){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  label = concentration = temperature = gene = value = colName = colNameOrig = 
    labelValue = foldChange = labelName = experiment <- NULL
  
  expInfo <- data %>% lapply(function(eSet){
    cbind(pData(eSet), colNameOrig = colnames(exprs(eSet)))
  }) %>%
    bind_rows(., .id = "experiment") %>%
    mutate(label = factor(as.character(label)))
  
  # Avoid problems if concentrations were imported from config file in 
  # scientific format. I.e. colname 'rel_fc_7.0000000000000007E-2' becomes 
  # 'rel_fc_7.0000000000000007E_2':
  expInfo <- expInfo %>% 
    mutate(colName = gsub("-", ".", expInfo$colNameOrig)) %>%
    mutate_if(is.character, factor)
  
  ## Assign concentration or temperature as label value
  if(any(grepl("concentration", colnames(expInfo)))){
    expInfo <- expInfo %>% rename(labelValue = concentration, labelName = label) 
  }
  if (any(grepl("temperature", colnames(expInfo)))){
    expInfo <- expInfo %>% rename(labelValue = temperature, labelName = label) 
  } 
  
  longExprs <- data %>% 
    purrr::map(biobroom::tidy.ExpressionSet) %>% 
    bind_rows(., .id = "experiment") %>%
    rename(id = gene, foldChange = value, colName = sample) %>%
    mutate_if(is.character, factor) %>%
    left_join(expInfo, by = c("experiment", "colName")) %>%
    select(-colName, colName = colNameOrig)
  
  ## If data was normalized, add suffix 'norm_' to the fold change column 
  ## names. Normalized data is recognized by the values of the normalization 
  ## coefficients in the fold change column annotation.
  flagIsNormalized <- any(!is.na(expInfo$normCoeff))
  
  if (flagIsNormalized) {
    longExprs <- longExprs %>% 
      mutate(colName = factor(paste("norm", colName, sep="_")))
  }
  
  ## Rearrange columns
  longExprs <- longExprs %>% 
    select(id, labelValue, foldChange, labelName, colName, experiment)
  
  return(longExprs)
}
