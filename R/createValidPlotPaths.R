createValidPlotPaths <- function(type, protIDs, plotDir){
  
  # type should be "meltCurve" or "smoothingSpline"
  
  ## File names for melting curve plots: Replace special characters by '_'.
  fNames <- paste0(type, "_", gsub("([^[:alnum:]])", "_", protIDs))
  
  # limit file name length to 255 characters to avoid crashes on most filesystems:
  maxLen <- 255 - nchar(".pdf") - nchar("_truncated") - nchar(as.character(length(fNames)))
  
  tooLong <- nchar(fNames) > maxLen
  
  cropSuffix <- paste0("_truncated", 1:sum(tooLong))
  
  fNames <- sapply(fNames, function(fTmp) {
    
    fNew <- substr(fTmp, 1, min(maxLen, nchar(fTmp)))
    
  }, simplify = TRUE, USE.NAMES = FALSE)
  
  fNames[tooLong] <- paste0(fNames[tooLong], cropSuffix)
  
  fNames <- paste0(fNames, ".pdf")
  
  plotPathsFull <- file.path(plotDir, fNames)
  
  pathTable <- data.frame(uniqueID = protIDs, path = plotPathsFull) %>%
    tbl_df() %>% mutate_all(as.character)
  
  return(pathTable)
}