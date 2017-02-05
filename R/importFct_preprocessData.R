importFct_preprocessData <- function(data, idVar){
  out <- data %>% 
    purrr::map(function(d) splitIDsIntoSeparateRows(singleDat = d, 
                                                    idVar = idVar))
  return(out)
}
