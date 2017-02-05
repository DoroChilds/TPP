extractConc <- function(configTable){
  # Return a vector of concentrations from a TPP-2D config table.
  # This vector will be used to create a TPP-CCR config table which enables
  # invoking 'analyzeTPPCCR' on the 2D dataset.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  label = concentration <- NULL
  
  allCols <- colnames(configTable)
  labelCols <- detectLabelColumnsInConfigTable(allColumns = allCols)
  
  uniqueConcentrations <- configTable %>% 
    subset(select = labelCols) %>%
    gather_("label", "concentration", labelCols) %>% 
    filter(concentration != "-") %>% 
    mutate(concentration = as.numeric(concentration)) %>% # Prevent sorting errors during CCR data import
    select(-label) %>%
    distinct %>% 
    extract2("concentration")
  
  return(uniqueConcentrations)
  
}
