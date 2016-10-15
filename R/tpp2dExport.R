#' @title Produce Excel table of 2D-TPP experiment.
#' @description Produce Excel table of 2D-TPP experiment analysis results.
#' 
#' 
#' @examples
#' data("panobinostat_2DTPP_smallExample")
#' load(system.file("example_data/2D_example_data/shortData2d.RData", package="TPP"))
#' # tpp2dExport(configTable = panobinostat_2DTPP_config, tab=shortData2d, 
#' #             resultPath=getwd(), 
#' #             idVar="representative", fcStr="norm_rel_fc_protein_", 
#' #             intensityStr="sumionarea_protein_", addCol=NULL)
#'
#' @return Creates excel file of the TPP-CCR analysis of the 2D-TPP data.
#' 
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param tab Table with results of the 2D-TPP analysis.
#' @param resultPath path for storing results table
#' @param idVar character string indicating how the column that contains the unique protein 
#'  identifiers is called
#' @param fcStr character string indicating how columns that contain the actual 
#'  fold change values are called
#' @param intensityStr character string indicating how columns that contain the raw masspec signal 
#'  (e.g. sumionareas) values are called
#' @param addCol additional names of columns which are to be attached to the result table
#' @param normalizedData boolean variabel indicating whether the data has been normalized
#' @param trRef character string containing a valid system path to a TPP-TR reference RData
#' file
#' 
#' 
#' @export
tpp2dExport <- function(configTable=NULL, tab=NULL, resultPath=NULL, idVar=NULL, fcStr=NULL, 
                        intensityStr=NULL, addCol=NULL, normalizedData=FALSE, trRef=NULL){
  if (is.null(fcStr)){
    fcStr <- "rel_fc_protein_"
  }
  fileName <- file.path(resultPath, paste(format(Sys.time(),"%Y-%m-%d"), "results_2D_TPP.xlsx", 
                                sep="_"))
  message("Writing results to file: ", fileName)
  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  boolPrefix <- c("passed_filter",
                  "protein_identified_in",
                  "sufficient_data_for_fit",
                  "model_converged",
                  "meets_FC_requirement", 
                  "pEC50_outside_conc_range")
  
  invisible(lapply(boolPrefix, function(bp){
    boolCol <- grep(bp, colnames(tab), value = TRUE)
    if (length(boolCol) > 0){
      newCol <- tab[,boolCol]
      newCol[which(newCol==TRUE)] <- "Yes"
      newCol[which(newCol==FALSE)] <- "No"
      newCol[which(is.na(newCol))] <- ""
      tab[,boolCol] <- newCol # new
    }
  }))
  
  ## Sort proteins in alphabetical order:
 # before continuing, make sure that the following command also works with dplyr::arrange (was written for plyr::arrange)
  tab <- arrange_(tab, .dots=c(idVar))
  
  # remove empty columns
  tab <- tab[, colSums(is.na(tab))<nrow(tab)]
  
  # find sumionarea and different fold change columns
  intensityCols <- colnames(tab)[which(grepl(intensityStr, colnames(tab)))]
  fcOrig <- colnames(tab)[which(grepl(paste("^", fcStr, sep=""), colnames(tab)))]
  fcTrans <- colnames(tab)[which(grepl(".*transformed", colnames(tab)))]
  lowestConc <- colnames(tab)[which(grepl(".*lowest_conc", colnames(tab)))]
  temperatureCol <- colnames(tab)[which(grepl("temperature", colnames(tab)))] # new
  
  # column oder
  if (normalizedData){
    fcNorm <- colnames(tab)[which(grepl(".*unmodified", colnames(tab)))]
    colOrder <- c(idVar, addCol, temperatureCol, intensityCols, fcOrig, fcNorm, fcTrans)
  }else{
    colOrder <- c(idVar, addCol, temperatureCol, intensityCols, fcOrig, fcTrans)
  }
  restCol <- setdiff(colnames(tab), c(colOrder, lowestConc))
  
  # rearrange columns of tab
  tab <- tab[, c(colOrder, restCol)]
  
  ## Generate plot-link columns for each protein
  tab <- tpp2dGeneratePlotLinkCols(tab, resultPath, idVar, trRef)
  allCols <- colnames(tab)
  
  # define excel styles
  headerS <- createStyle(fgFill = "#DCE6F1", border="Bottom", textDecoration="Bold")
  darkgreen <- createStyle(bgFill="darkgreen")
  lightgreen <- createStyle(bgFill="darkolivegreen3")
  khaki <- createStyle(bgFill="khaki1")
  lightorange <- createStyle(bgFill="goldenrod1")
  darkorange <- createStyle(bgFill="darkorange")
  fc_cols <- grep(".*unmodified", colnames(tab))
  
  # check whether any colnames are duplicated
  if(any(duplicated(tolower(colnames(tab))))){
    colnames(tab)[which(duplicated(tolower(colnames(tab))))] <- 
      paste0(colnames(tab)[which(duplicated(tolower(colnames(tab))))], "_2")
  }
  
  # create workbook and fill with table columns:
  wb <- createWorkbook()
  addWorksheet(wb, "Exp details")
  writeDataTable(wb, sheet="Exp details",  x=configTable, startCol=1, startRow=1, 
                 rowNames=FALSE, colNames=TRUE, headerStyle=headerS)
  addWorksheet(wb, "pEC50")
  writeDataTable(wb, sheet="pEC50",  x=tab, startCol=1, startRow=1, 
                 rowNames=FALSE, colNames=TRUE, headerStyle=headerS)
  conditionalFormatting(wb, sheet="pEC50", cols=fc_cols, rows=2:(nrow(tab)+1), rule=">=1.5", style=lightgreen)
  conditionalFormatting(wb, sheet="pEC50", cols=fc_cols, rows=2:(nrow(tab)+1), rule="<=0.5", style=darkorange)
  conditionalFormatting(wb, sheet="pEC50", cols=fc_cols, rows=2:(nrow(tab)+1), rule="<=0.67", style=lightorange)
  conditionalFormatting(wb, sheet="pEC50", cols=fc_cols, rows=2:(nrow(tab)+1), rule="<1.5", style=khaki)
  conditionalFormatting(wb, sheet="pEC50", cols=fc_cols, rows=2:(nrow(tab)+1), rule=">=2", style=darkgreen) # last rule appears in Excel as first rule
  
  
  ## Add column with links to fitted curve plots:
  linkCol <- grepl("plot_", allCols)
  if (any(linkCol)){
    invisible(sapply(colnames(tab)[which(linkCol)], function(pp){
      ind <- which(colnames(tab)==pp)
      class(tab[,pp]) <- 'hyperlink'
      suppressWarnings(writeData(wb, sheet="pEC50", x=tab[,pp], startCol=ind, 
                               startRow = 2, keepNA=FALSE))
    }))
  }
  
  tryCatch({
    saveWorkbook(wb, file=fileName, overwrite=TRUE)
    message("File created successfully!\n")
  },
  error = function(err){
    message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
    message(err)
    message(paste("\n\nAlthough the Excel output failed, you can still access the results of function 
                  'tpptrAnalyzeMeltingCurves' via its return value or as an R object in the results 
                  folder.\n"))
  })
}
