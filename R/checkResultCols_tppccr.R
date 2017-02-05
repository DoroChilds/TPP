checkResultCols_tppccr <- function(curveParDF, fcFilterDF, minR2){
  ## Retrieve relevant columns for quality check:
  r2Cols <- grep("R_sq", colnames(curveParDF), value=TRUE)
  r2     <- subset(curveParDF, select = r2Cols)
  
  fcReqCols <- grep("meets_FC_requirement", colnames(fcFilterDF), value=TRUE)
  fcReq     <- subset(fcFilterDF, select = fcReqCols)
  
  ##  Check which protein passed the filter criteria
  passed_filter_r2 <-r2>=minR2
  passed_filter_fc <-fcReq
  
  passed_filter_both <- passed_filter_r2 & passed_filter_fc
  # Replace NAs by FALSE. NAs can appear, for example, in proteins without
  # melting curve fits, and hence no R2 values.
  passed_filter_both[is.na(passed_filter_both)] <- FALSE
  
  newCols <- gsub("R_sq", replacement="passed_filter", r2Cols)
  colnames(passed_filter_both) <- newCols
  
  outDF <-data.frame("Protein_ID"=curveParDF$Protein_ID, passed_filter_both)
  return(outDF)
}