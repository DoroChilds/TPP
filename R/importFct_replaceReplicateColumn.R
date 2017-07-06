importFct_replaceReplicateColumn <- function(cfg){
  ## This function checks whether table contains the deprecated column 
  ## 'Replicate', and assigns valid comparisons instead.
  ## In the first version of the TPP package, replicate columns were assigned
  ## to derive the comparisons of interest (vehicle replicates 1 & 2 vs. 
  ## treatment replicates 1 & 2). Later, they were replaced by comparison
  ## columns to increase flexibility in defining the contrast use for testing.
  ## This function ensures backwards compatibility by converting an `old-style'
  ## config table as described in the original publication, to the modern 
  ## version.
  
  
  ## Look for replicate column:
  replicate_column <- colnames(cfg)[colnames(cfg) =="Replicate"]
  comparison_columns <- grep("Comparison", colnames(cfg), value = TRUE)
  
  if (length(replicate_column) > 0) {
    
    ## Assign comparison columns:
    if (length(comparison_columns) == 0){
      
      message("Replacing config table column 'Replicate' by corresponding columns starting with suffix 'Comparison'.\n")
      
      replicates <- cfg[[replicate_column]]
      rep_levels <- unique(replicates)
      
      for (r in rep_levels) {
        new_col <- ifelse(replicates == r, "x", "")
        new_col_name <- paste0("Comparison", "VT", r)
        cfg[[new_col_name]] <- new_col
      }
      
      ## Remove obsolete replicate column:
      cfg[[replicate_column]] <- NULL
      
    } else {
      
      ## If both column types were present, throw an informative error message
      ## to avoid unexpected effects in the downstream test results:
      stop("Cannot assign contrasts from config table. Config table may have either a column named 'Replicate', or columns starting with prefix 'Comparison'but not both!")
      
    }
  } 
  
  return(cfg)
}
