convert_2dData_wide_to_long <- function(datWide, idColname, fcStr){
  # Convert 2D-TPP dataset to long table
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  columnName = fc = uniqueID <- NULL
  
  # ptrn <- paste0(idColname, 
  #                "|temperature|", 
  #                "^", fcStr, "[0-9,\\.]+$|", 
  #                "^", fcStr, "[0-9,\\.]+_unmodified")
  ptrn <- paste0(idColname,"|temperature|", fcStr)
  datLong <- datWide %>% tibble::as_tibble() %>%
    select(matches(ptrn)) %>% 
    gather(columnName, fc, contains(fcStr)) %>%
    rename_(uniqueID = idColname) %>%
    arrange(uniqueID)
  
  # Add column with drug concentrations
  ptrn <- paste(fcStr, "([0-9,\\.]+)[^0-9]*", sep="")
  oldValues <- unique(datLong$columnName)
  newValues <- paste(sub(ptrn, "\\1", oldValues), "uM", sep="")
  newCol <- plyr::mapvalues(datLong$columnName, oldValues, newValues)
  datLong2 <- datLong %>% mutate(drugConc = newCol)
}