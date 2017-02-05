exportFct_ensureUniqueColumns <- function(dat){
  # check whether any colnames are duplicated
  chkCols <- colnames(dat) %>% tolower
  isDupl <- duplicated(chkCols)
  if(any(isDupl)){
    newCols <- ifelse(isDupl, paste0(chkCols, "_2"), chkCols)
    colnames(dat) <- newCols
  }
  return(dat)
}
