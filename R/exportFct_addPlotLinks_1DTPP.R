exportFct_addPlotLinks_1DTPP <- function(wb, sheet, dat){
  ## Add column with links to fitted curve plots
  ## to do: combine this function with 'exportFct_addPlotLinks_2DTPP'
  allCols <- colnames(dat)
  
  for (plotCol in c("splinefit_plot", "meltcurve_plot")){
    linkCol <- grepl(plotCol, allCols)
    if (any(linkCol)){
      relPaths <- as.character(dat[,linkCol])
      plotPaths <- relPaths
      plotPaths <- ifelse(is.na(plotPaths), "none", file.path(".", plotPaths))
      names(plotPaths) <- ifelse(plotPaths=="none", "", 
                                 gsub("([^[:alnum:]])", "_", dat$Protein_ID))
      class(plotPaths) <- "hyperlink"
      #plotPaths[plotPaths=="_"] <- NA
      suppressWarnings(writeData(wb, sheet=sheet, x=plotPaths, startCol=which(linkCol), 
                                 startRow = 2, keepNA = TRUE))
    }    
  }

  return(wb)
}
