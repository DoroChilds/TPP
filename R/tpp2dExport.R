#' @title Produce Excel table of 2D-TPP experiment.
#' @description Produce Excel table of 2D-TPP experiment analysis results.
#' 
#' 
#' @examples
#' data(panobinostat_2DTPP_smallExample)
#' load(system.file("example_data/2D_example_data/shortData2d.RData", package="TPP"))
#' # tpp2dExport(configTable = panobinostat_2DTPP_config, tab=shortData2d, 
#' #             outPath=getwd(), 
#' #             idVar="representative", fcStr="norm_rel_fc_protein_", 
#' #             intensityStr="sumionarea_protein_", addCol=NULL)
#' 
#' data(panobinostat_2DTPP_smallExample)
#' # cfgRaw <- panobinostat_2DTPP_config
#' # datRaw <- panobinostat_2DTPP_data
#' # datIn <- tpp2dImport(cfgIn, datRaw, fcStr = NULL)
#' # datFC <- tpp2dComputeFoldChanges(data = datIn)
#' # datNorm <- tpp2dNormalize(data = datFC)
#' # cfgCCR <- convert_2D_cfgTable_to_CCR_cfgTable(cfgIn)
#' # datFitted <- tpp2dCurveFit(datNorm, nCores = 2)
#' # tpp2dCreateReport(getwd(), cfgIn, resultTable = datFitted, idVar = "representative", 
#' #                   intensityStr = "sumionarea_protein_")
#' # tpp2dExport(tab = datFitted, outPath = getwd(), addPlotColumns = FALSE) 
#'
#' @return Creates excel file of the TPP-CCR analysis of the 2D-TPP data.
#' 
#' @param tab Table with results of the 2D-TPP analysis.
#' @param outPath path for storing results table
#' @param addCol additional names of columns which are to be attached to the result table
#' @param trRef character string containing a valid system path to a TPP-TR reference RData
#' file
#' @param addPlotColumns boolean variable indicating whether paths to plot files 
#' should be generated and checked for validity. De-activate if no dose-response
#' curve plots were produced during the analysis.
#' @param configTable DEPRECATED
#' @param resultPath DEPRECATED
#' @param idVar DEPRECATED
#' @param fcStr DEPRECATED
#' @param intensityStr DEPRECATED
#' @param normalizedData DEPRECATED
#' 
#' 
#' @export
tpp2dExport <- function(configTable = NULL, 
                        tab, 
                        resultPath = NULL, 
                        idVar = NULL, 
                        fcStr = NULL, 
                        intensityStr = NULL, 
                        outPath, 
                        addCol = NULL, 
                        normalizedData = NULL, 
                        trRef = NULL, 
                        addPlotColumns = TRUE){
  
  if (!missing(configTable)){
    warning("`configTable` is deprecated.", call. = TRUE)
  }
  
  if (!missing(resultPath)){
    warning("`resultPath` is deprecated.", call. = TRUE)
  }
  
  if (!missing(idVar)){
    warning("`idVar` is deprecated.", call. = TRUE)
  }
  
  if (!missing(fcStr)){
    warning("`fcStr` is deprecated.", call. = TRUE)
  }
  
  if (!missing(intensityStr)){
    warning("`intensityStr` is deprecated.", call. = TRUE)
  }
  
  if (!missing(normalizedData)){
    warning("`normalizedData` is deprecated.", call. = TRUE)
  }
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("tab", "outPath"))
  
  # Obtain config table used for data import (stored as attribute of imported data):
  configTable <- attr(tab, "configTable")
  
  # Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(tab, "importSettings")
  idVar <- checkAndReturnDataSetting(importSettings, "proteinIdCol", colnames(tab))
  fcStr <- checkAndReturnDataSetting(importSettings, "fcStr", colnames(tab))
  intensityStr <- checkAndReturnDataSetting(importSettings, "intensityStr", colnames(tab))
  normalizedData <- !is.null(importSettings$fcStrNorm)
  
  
  fTmp <- paste0(format(Sys.time(),"%Y-%m-%d"), "_results_2D_TPP.xlsx")
  fileName <- file.path(outPath, fTmp)
  message("Writing results to file: ", fileName)
  
  # Boolean columns: Convert TRUE/FALSE to "Yes"/"No" values
  tab <- exportFct_convertBoolean_2DTPP(tab)
  
  ## Sort proteins in alphabetical order:
  tab <- arrange_(tab, .dots=c(idVar))
  
  # remove empty columns
  # tab <- tab[, colSums(is.na(tab))<nrow(tab)]
  
  ## Rearrange columns:
  tab <- exportFct_sortCols(dat = tab, 
                            idVar = idVar, addCol = addCol, 
                            intensityStr = intensityStr, fcStr = fcStr, 
                            normalizedData = normalizedData)
  
  ## Generate plot-link columns for each protein
  if (addPlotColumns){
    tab <- exportFct_addPlotColumns(tab = tab, path = outPath, idVar = idVar, 
                                    trRef = trRef)
  }
  
  # check whether any colnames are duplicated
  tab <- exportFct_ensureUniqueColumns(tab)
  
  ## Convert decimal separators into the format used by the current OS in the 
  ## label columns. This is currently necessary because they are of 
  ## class 'character' due to the '-' entries which are present by default
  ## in 2D-config tables.
  allCfgCols     <- colnames(configTable)
  
  labelCols <- detectLabelColumnsInConfigTable(allColumns = allCfgCols)
  
  configTable[,labelCols] <- suppressWarnings(
    apply(configTable[,labelCols], 2, as.numeric)
  )
  
  # create workbook and store data:
  wb <- createWorkbook()
  addWorksheet(wb, "Exp details")
  addWorksheet(wb, "pEC50")
  headerStyle <- createStyle(fgFill = "#DCE6F1", border="Bottom", textDecoration="Bold")
  
  writeDataTable(wb, sheet = "Exp details", x = configTable, 
                 startCol = 1, startRow = 1, 
                 rowNames = FALSE, colNames = TRUE, 
                 headerStyle = headerStyle)
  
  writeDataTable(wb, sheet = "pEC50",  x = tab, 
                 startCol = 1, startRow = 1, 
                 rowNames = FALSE, colNames = TRUE, 
                 headerStyle = headerStyle)
  
  ## Color-code fold changes:
  #  if (normalizedData){
  fc_cols <- grep(".*unmodified", colnames(tab))
  #  } else {
  #    fc_cols <- grep(paste("^", fcStr, sep=""), allCols)
  #  }
  wb <- exportFct_colorCodeColumns(wb = wb, sheet = "pEC50", 
                                   cols = fc_cols, dat = tab)
  
  ## Convert plot column entries to links:
  wb <- exportFct_addPlotLinks_2DTPP(wb = wb, sheet = "pEC50", dat = tab)
  
  ## Save final table
  success <- exportFct_trySave(wb = wb, file = fileName)
  return(fileName)
}
