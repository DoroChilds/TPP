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
  message("Writing results to file: ", file)
  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  allCols <- colnames(tab)
  boolPrefix <- c("passed_filter",
                  "protein_identified_in",
                  "sufficient_data_for_fit",
                  "model_converged",
                  "min_pVals_less_0.05_and_max_pVals_less_0.1",
                  "meltP_diffs_have_same_sign",
                  "meltP_diffs_T_vs_V_greater_V1_vs_V2",
                  "minSlopes_less_than_0.06",
                  "fulfills_all_4_requirements",
                  "meets_FC_requirement", 
                  "pEC50_outside_conc_range")
  
  for (bp in boolPrefix){
    boolCols <- grep(bp, colnames(tab), value = TRUE)
    for (bc in boolCols){
      x <-tab[,bc]
      xNew <- rep(NA_character_, length(x))
      xNew[which(x==TRUE)] <- "Yes"
      xNew[which(x==FALSE)] <- "No"
      xNew[which(is.na(x))] <- ""
      tab[,bc] <- xNew      
    }
  }
  
  ## Remove inflection point column from TPP-TR output:
  rmCols <- c(grep("inflPoint", allCols, value=TRUE),
              allCols[substr(allCols, 1, 2) == "a_"],
              allCols[substr(allCols, 1, 2) == "b_"])
  tab <- subset(tab, select=!(allCols %in% rmCols))
  allCols <- colnames(tab)
  
  ## Remove plot column from TPP-TR output if it only contains missing values:
  if (all(is.na(tab$plot))){
    tab <- subset(tab, select=!(allCols %in% "plot"))    
    allCols <- colnames(tab)
  }
  
  
  ## Sort proteins in alphabetical order:
  tab <- arrange(df=tab, tab$Protein_ID)
  
  ## Create workbook and fill with table columns:
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  writeDataTable(wb, sheet="Results",  x=tab, startCol=1, startRow=1, 
                 rowNames=FALSE, colNames=TRUE)
  
  ## Add column with links to fitted curve plots:
  linkCol <- grepl("plot", allCols)
  if (any(linkCol)){
    relPaths <- as.character(tab[,linkCol])
    plotPaths <- relPaths
    plotPaths <- ifelse(is.na(plotPaths), "none", file.path(".", plotPaths))
    names(plotPaths) <- ifelse(plotPaths=="none", "", gsub("([^[:alnum:]])", "_", tab$Protein_ID))
    class(plotPaths) <- "hyperlink"
    #plotPaths[plotPaths=="_"] <- NA
    suppressWarnings(writeData(wb, sheet="Results", x=plotPaths, startCol=which(linkCol), 
                               startRow = 2, keepNA=TRUE))
  }
  
  
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
  
  tryCatch({
    saveWorkbook(wb, file=file, overwrite=TRUE)
    message("File created successfully!\n")
  },
  error = function(err){
    message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
    message(err)
    message(paste("\n\nAlthough the Excel output failed, you can still access the results of function 'tpptrAnalyzeMeltingCurves' via its return value or as an R object in the results folder.\n"))
  })
}
