exportFct_addPlotLinks_2DTPP <- function(wb, sheet, dat){
  ## Add column with links to fitted curve plots
  ## to do: combine this function with 'exportFct_addPlotLinks_1DTPP'
  allCols <- colnames(dat)
  ## Add column with links to fitted curve plots:
  linkCol <- grepl("plot_", allCols)
  if (any(linkCol)){
    invisible(sapply(colnames(dat)[which(linkCol)], function(pp){
      ind <- which(allCols == pp)
      class(dat[,pp]) <- 'hyperlink'
      suppressWarnings(writeData(wb, sheet=sheet, x=dat[,pp], startCol=ind, 
                                 startRow = 2, keepNA = FALSE))
    }))
  }
  return(wb)
}
