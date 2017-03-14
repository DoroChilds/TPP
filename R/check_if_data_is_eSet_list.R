check_if_data_is_eSet_list <- function(data){
  
  data_is_list <- is.list(data) & !is.data.frame(data)
  
  if (data_is_list){
    data_is_eSet_list <- sapply(data, class) %>% unique %>% 
      identical("ExpressionSet")
  } else {
    data_is_eSet_list <- FALSE
  }
  
  return(data_is_eSet_list)
  
}