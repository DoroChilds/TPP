determine_data_format <- function(data){
  
  data_is_eSet_list <- check_if_data_is_eSet_list(data)
  
  data_is_dataframe <- is.data.frame(data)
  
  if (data_is_dataframe){
    
    data_is_wide_table <- data[[idColumn]] %>% 
      table %>% 
      unique %>% 
      identical(1)
    
    data_is_long_table <- !data_is_wide_table
    
  } else {
    data_is_wide_table <- FALSE
    data_is_long_table <- FALSE
  }
  
  if (data_is_eSet_list){
    out <- "eSetList"
  } else if (data_is_long_table){
    out <- "longTable"
  } else if (data_is_wide_table){
    out <- "wideTable"
  } 
  
  return(out)
}
