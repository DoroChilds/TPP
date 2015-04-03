#' @title Produce Excel table of TPP-TR or TPP-CCR experiment.
#' 
#' @examples
#' data(hdacTR_resultsTable_smallExample)
#' tppExport(resultTable, "tpptr_example_results.xlsx")
#' 
#' @param tab Table with results of the TPP analysis.
#' @param file path for storing results table
#' @return No value returned.
#' @export
tppExport <- function(tab, file){
  message("Writing results to file: ", file)
  # Boolean columns: Convert TRUE/(FALSE|NA) to YES/NO values
  allCols <- colnames(tab)
  boolCols <- c(grep("model_converged_", allCols, value=TRUE),
                grep("sufficient_data_for_fit_", allCols, value=TRUE),
                grep("passed_filter", allCols, value=TRUE),
                grep("min_pVals_less_0.05_and_max_pVals_less_0.1", allCols, 
                     value=TRUE),
                grep("meltP_diffs_have_same_sign", allCols, value=TRUE),
                grep("meltP_diffs_T_vs_V_greater_V1_vs_V2", allCols, value=TRUE),
                grep("minSlopes_less_than_0.06", allCols, value=TRUE),
                grep("fulfills_all_4_requirements", allCols, value=TRUE))
  for (bc in boolCols){
    x <-tab[,bc]
    xNew <- rep(NA_character_, length(x))
    xNew[x==TRUE] <- "Yes"
    xNew[x==FALSE | is.na(x)] <- "No"
    tab[,bc] <- xNew
  }
  
  ## Remove inflection point column from TPP-TR output:
  inflPointColumn <- grep("inflPoint", allCols, value=TRUE)
  tab <- subset(tab, select=!(allCols %in% inflPointColumn))
  allCols <- colnames(tab)
  
  ## Sort proteins in alphabetical order:
  tab <- arrange(df=tab, tab$Protein_ID)
  
  ## Create workbook and fill with table columns:
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  writeDataTable(wb, sheet="Results",  x=tab, startCol=1, startRow=1, 
                 rowNames=FALSE, colNames=TRUE)
  
  ## Add column with links to fitted curve plots:
  linkCol <- grep("Plot", allCols)
  plotLinks <- as.character(tab[,linkCol])
  names(plotLinks) <- ifelse(is.na(plotLinks), "", as.character(tab$Protein_ID))
  plotLinks[is.na(plotLinks)] <- "none"  
  class(plotLinks) <- "hyperlink"
  suppressWarnings(writeData(wb, sheet="Results", x=plotLinks, startCol=linkCol, 
                             startRow = 2))
  
  tryCatch({
    saveWorkbook(wb, file=file, overwrite=TRUE)
    message("File created successfully!\n")
  },
  error = function(err){
    message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
    message(err)
    message(paste("\n\nAlthough the Excel output failed, you can still access the results of function 'tpptrResultTable' via its return value or as an R object in the results folder.\n"))
  })
}