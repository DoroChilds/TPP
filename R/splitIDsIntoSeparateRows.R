splitIDsIntoSeparateRows <- function(singleDat, idVar) {
  # cat(paste(Sys.time(), "checking for samesetproteins that also are representatives...\n"))
  ids <- singleDat[[idVar]]
  
  singleDat$row_derived_from_non_unique_identifiers <- FALSE
  
  nonUniques <- grep("\\|", ids, value = TRUE)
  
  newRows <- lapply(nonUniques, function(entry){
    ## Create one row for each unique identifier in 'entry'
    currentRow <- singleDat %>% filter_(paste0(idVar, '  == "', entry, '"'))
    
    # Repeat the row corresponding to the current entry separately for each unique
    # identifier in the entry.  
    uniques <- unique(unlist(strsplit(entry, "\\|"))) %>% sort
    newRows <-  currentRow %>% colwise(function(x) rep(x, length(uniques)))(.)
    newRows[[idVar]] <- uniques
    return(newRows)
  }) %>%
    bind_rows %>% 
    mutate(row_derived_from_non_unique_identifiers = TRUE)
  
  # # check for each protein in a split column if it already occurs in the data set
  # new_rows <- data.frame()
  # uniqueIDs <- unique(unlist(strsplit(ids, "\\|")))
  # for (p in uniqueIDs) {
  #   l_index <- grep(p, ids)
  #   l_num <- length(l_index)
  #   if (l_num > 0) {
  #     cat(paste(p, l_num, "\n"))
  #     new <- singleDat[l_index,]
  #     new[[idVar]] <- p
  #     # new$clustername <- a_results[grep(p, a_results$representative)[1],"clustername"]
  #     
  #     # add new rows
  #     new_rows <- rbind(new_rows, new)
  #   }
  # }
  
  
  if (nrow(newRows) > 0) {
    message("Splitting ", length(nonUniques), 
            " rows with multiple identifiers into ", nrow(newRows), 
            " separate rows.")
  } 
  newDat <- singleDat[!singleDat[[idVar]] %in% nonUniques, ] %>%
    rbind(newRows)
  
  
  return(newDat)
}
