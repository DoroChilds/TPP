detectLabelColumnsInConfigTable <- function(allColumns){
  ## Find and return the column names belonging to isobaric labels in a 
  ## user-defined config table in wide format.
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("allColumns"))
  
  ## Comparison columns (used for TR-config tables):
  compCols    <- grep("comparison", allColumns, value=TRUE, ignore.case=TRUE)
  
  ## Find all columns NOT belonging to an isobaric label:
  noLabelCols <- c("Experiment", "Path", # General columns
                   "Condition", compCols, # TR-specific columns
                   "Compound", "Temperature", "RefCol") # 2D-TPP specific columns
  
  ## Detect the columns belonging to isobaric labels:
  labelCols <- setdiff(allColumns, noLabelCols)
  
  return(labelCols)
}