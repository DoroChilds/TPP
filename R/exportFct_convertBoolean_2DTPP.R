exportFct_convertBoolean_2DTPP <- function(dat){
  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  # to do: combine this function with 'exportFct_convertBoolean_1DTPP'
  boolPrefix <- c("passed_filter",
                  "protein_identified_in",
                  "sufficient_data_for_fit",
                  "model_converged",
                  "meets_FC_requirement", 
                  "pEC50_outside_conc_range")
  
  for (bp in boolPrefix){
    boolCols <- grep(bp, colnames(dat), value = TRUE)
    for (bc in boolCols){
      x <- dat[,bc]
      xNew <- rep(NA_character_, length(x))
      xNew[which(x==TRUE)] <- "Yes"
      xNew[which(x==FALSE)] <- "No"
      xNew[which(is.na(x))] <- ""
      dat[,bc] <- xNew      
    }
  } 
  return(dat)
}
