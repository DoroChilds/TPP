obtain_fcStr_from_df_annotation <- function(dat){
  message("Checking which columns in the data table contain the fold change values for fitting and plotting...")
  importSettings <- attr(dat, "importSettings")
  
  if (is.null(importSettings)){
    msg0 <- paste0("Data table needs an attribute named 'importSettings'. ",
                   "This attribute is a list containing all settings that ",
                   "were used when importing and processing the data.")
    stop(msg0)
  }
  
  fcStr     <- importSettings$fcStr
  fcStrNorm <- importSettings$fcStrNorm
  
  ## Assess whether the data was normalized (the field 'fcStrNorm' was 
  ## filled by the normalization function and is empty otherwise):
  if (!is.null(fcStrNorm)){
    finalFcPrefix <- fcStrNorm
    msg1 <- paste0("Normalized data columns detected with prefix '", 
                   finalFcPrefix,"'. Analysis will be based on these values.")
  } else {
    finalFcPrefix <- fcStr
    msg1 <- paste0("No normalized data columns detected. Analysis will be 
                   based on the unnormalized values obtained from columns 
                   with prefix '", finalFcPrefix, "'.") %>%
      gsub("\n *", "", .)
  }
  msg2 <- paste0("This information was found in the attributes of the input data 
                 (access with attr(dataTable, 'importSettings'))") %>%
    gsub("\n *", "", .)
  
  ## Test if the chosen prefix really exists in the table columns:
  allCols <- colnames(dat)
  validStr <- grepl(finalFcPrefix, allCols) %>% any
  if (!validStr){
    stop("The detected fold change column prefix '", finalFcPrefix, 
         "' could not be found in the columns of the data table.")
  } else {
    message(msg1, "\n", msg2)
    return(finalFcPrefix)
  }
}
