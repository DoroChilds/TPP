exportFct_convertBoolean_1DTPP <- function(dat){
  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  # to do: combine this function with 'exportFct_convertBoolean_2DTPP'
  boolPrefix <- c("passed_filter",
                  "protein_identified_in",
                  "sufficient_data_for_fit",
                  "model_converged",
                  "min_pVals_less_0.1_and_max_pVals_less_0.2",
                  "meltP_diffs_have_same_sign",
                  "meltP_diffs_T_vs_V_greater_V1_vs_V2",
                  "minSlopes_less_than_0.06",
                  "fulfills_all_4_requirements",
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
