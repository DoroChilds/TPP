#' @title Produce Excel table of TPP-TR or TPP-CCR experiment.
#' @description Produce Excel table of TPP-TR or TPP-CCR experiment out of the
#' data frame returned by \code{\link{tpptrAnalyzeMeltingCurves}}
#' 
#' @examples
#' data(hdacTR_resultsTable_smallExample)
#' tppExport(resultTable, "tpptr_example_results.xlsx")
#' 
#' @param tab Table with results of the TPP analysis.
#' @param file path for storing results table
#' @param expColors character vector of background colors to group the result 
#' columns belonging to different experiments.
#' @param expNames character vector of experiment names of the same length as 
#' expColors. 
#' @return No value returned.
#' @export
tppExport <- function(tab, file, expNames=NULL, expColors=NULL){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Protein_ID <- NULL
  
  message("Writing results to file: ", file)
  allCols <- colnames(tab)

  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  tab <- exportFct_convertBoolean_1DTPP(tab)
  
  ## Remove inflection point column from TPP-TR output:
  rmCols <- c(grep("inflPoint", allCols, value=TRUE),
              allCols[substr(allCols, 1, 2) == "a_"],
              allCols[substr(allCols, 1, 2) == "b_"])
  tab <- subset(tab, select=!(allCols %in% rmCols))
  allCols <- colnames(tab)
  
  ## Remove plot column from TPP-TR output if it only contains missing values:
  for (plotCol in c("splinefit_plot", "meltcurve_plot")){
    if (!any(grepl(plotCol, colnames(tab)))){
      tab[[plotCol]] <- NA
    }
    if (all(is.na(tab[[plotCol]]))){
      tab[[plotCol]] <- NULL
    } 
  }
  allCols <- colnames(tab)
  
  ## Sort proteins in alphabetical order:
  tab <- arrange(tab, Protein_ID)
  
  ## Create workbook and fill with table columns:
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  writeDataTable(wb, sheet="Results",  x=tab, startCol=1, startRow=1, 
                 rowNames=FALSE, colNames=TRUE)
  
  ## Add column with links to fitted curve plots:
  wb <- exportFct_addPlotLinks_1DTPP(wb = wb, sheet = "Results", dat = tab)
  
  
  ## Mark experiments by different colors in the result table:
  if (!is.null(expNames) && !is.null(expColors)){
    if (length(expNames) == length(expColors)){
      for (i in 1:length(expNames)){
        cTmp <- expColors[i]
        nTmp <- expNames[i]
        j <-grep(nTmp, colnames(tab))
        s <- createStyle(fgFill=cTmp)
        addStyle(wb=wb, sheet=1, style = s, rows=1:(nrow(tab)+1), cols=j, gridExpand = TRUE, stack=TRUE)
      }
    } else {
      warning("Arguments 'expNames' and 'expColors' need to have same length. Cannot group columns belonging to the same experiment by color.")
    }
  }
  
  success <- exportFct_trySave(wb = wb, file = file)
  return(success)
}
