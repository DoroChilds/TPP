detectLabelColumnsInConfigTable <- function(allColumns){
  ## Find and return the column names belonging to isobaric labels in a 
  ## user-defined config table in wide format.
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("allColumns"))
  
  ## Find all columns NOT belonging to an isobaric label:
  noLabelCols <- nonLabelColumns()$column %>% as.character %>% unique
  
  ## Comparison columns (used for TR-config tables):
  # to do: ignore case for every column name to ensure consistency in treatment 
  # of the different columns
  compCols    <- grep("comparison", allColumns, value=TRUE, ignore.case=TRUE)
  noLabelCols <- c(noLabelCols, compCols)
  
  ## Detect the columns belonging to isobaric labels:
  labelCols <- setdiff(allColumns, noLabelCols)
  
  return(labelCols)
}

nonLabelColumns <- function(){
  
  out <- data.frame(
    column = c("Experiment", "Experiment", "Experiment",
               "Path", "Path", "Path", # General columns
               "Condition", "Replicate", # TR-specific columns
               "Compound", "Temperature", "RefCol"), # 2D-TPP specific columns
    type = c("TR", "CCR", "2D",
             "TR", "CCR", "2D",
             "TR", "TR", 
             "2D", "2D", "2D"),
    obligatory = c(TRUE, TRUE, TRUE, 
                   FALSE, FALSE, FALSE, 
                   TRUE, FALSE, 
                   TRUE, TRUE, TRUE),
    exclusive = c(FALSE, FALSE, FALSE,
                  FALSE, FALSE, FALSE,
                  TRUE, TRUE,
                  TRUE, TRUE, TRUE)
  )
  
  return(out)
}