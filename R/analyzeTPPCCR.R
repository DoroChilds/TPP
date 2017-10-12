#' @title Analyze TPP-CCR experiment
#' @description Performs analysis of a TPP-CCR experiment by invoking routines 
#'   for data import, data processing, normalization, curve fitting, and 
#'   production of the result table.
#'   
#' @return A data frame in which the fit results are stored row-wise for each
#'   protein.
#'   
#' @references Savitski, M. M., Reinhard, F. B., Franken, H., Werner, T., 
#'   Savitski, M. F., Eberhard, D., ... & Drewes, G. (2014). Tracking cancer 
#'   drugs in living cells by thermal profiling of the proteome. Science, 
#'   346(6205), 1255784.
#'   
#'   Franken, H, Mathieson, T, Childs, D. Sweetman, G. Werner, T. Huber, W. & Savitski, M. M. (2015),
#'   Thermal proteome profiling for unbiased identification of drug targets and detection of downstream effectors.
#'   Nature protocols 10(10), 1567-1593.
#'   
#' @details Invokes the following steps: \enumerate{ \item Import data using the
#'   \code{\link{tppccrImport}} function. \item Perform normalization by fold 
#'   change medians (optional) using the \code{\link{tppccrNormalize}} function.
#'   To perform normalization, set argument \code{normalize=TRUE}. \item Fit and
#'   analyze dose response curves using the \code{\link{tppccrCurveFit}} 
#'   function. \item Export results to Excel using the \code{\link{tppExport}} 
#'   function. }
#'   
#'   The default settings are tailored towards the output of the python package 
#'   isobarQuant, but can be customized to your own dataset by the arguments 
#'   \code{idVar, fcStr, naStrs, qualColName}.
#'   
#'   If \code{resultPath} is not specified, result files are stored at the path 
#'   defined in the first entry of \code{configTable$Path}. If the input data are not 
#'   specified in \code{configTable}, no result path will be set. This means 
#'   that no output files or dose response curve plots are produced and 
#'   \code{analyzeTPPCCR} just returns the results as a data frame.
#'   
#'   The function \code{analyzeTPPCCR} reports intermediate results to the 
#'   command line. To suppress this, use \code{\link{suppressMessages}}.
#'   
#' The dose response curve plots will be stored in a subfolder with 
#'   name \code{DoseResponse_Curves} at the location specified by 
#'   \code{resultPath}.
#'   
#' Only proteins with fold changes bigger than
#' \code{[fcCutoff * (1 - fcTolerance)} or smaller than 
#' \code{1/(fcCutoff * (1 - fcTolerance))]} will be used for curve fitting.
#' Additionally, the proteins fulfilling the fcCutoff criterion without 
#' tolerance will be marked in the output column \code{meets_FC_requirement}.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrResults <- analyzeTPPCCR(configTable=hdacCCR_config, 
#'                                data=hdacCCR_data, nCores=1)
#'   
#' @param configTable dataframe, or character object with the path to a file, 
#'   that specifies important details of the TPP-CCR experiment. See Section 
#'   \code{details} for instructions how to create this object.
#' @param data single dataframe, containing fold change measurements and 
#'   additional annotation columns to be imported. Can be used instead of 
#'   specifying the file path in the \code{configTable} argument.
#' @param resultPath location where to store dose-response curve plots and 
#'   results table.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param normalize perform median normalization (default: TRUE).
#' @param ggplotTheme ggplot theme for dose response curve plots.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param plotCurves boolean value indicating whether dose response curves should
#'   be plotted. Deactivating plotting decreases runtime.
#' @param verbose print name of each fitted or plotted protein to the command 
#' line as a means of progress report.
#' @param xlsxExport produce results table in xlsx format and store at the 
#' location specified by the \code{resultPath} argument.
#'   
#' @seealso tppDefaultTheme
#'   
#' @export
analyzeTPPCCR <- function(configTable, data=NULL, resultPath=NULL, 
                          idVar="gene_name", fcStr="rel_fc_", 
                          naStrs=c("NA", "n/d", "NaN", "<NA>"), 
                          qualColName="qupm",
                          normalize=TRUE, ggplotTheme=tppDefaultTheme(), 
                          nCores="max", nonZeroCols="qssm", 
                          r2Cutoff=0.8,  fcCutoff=1.5, slopeBounds=c(1,50),
                          plotCurves=TRUE, verbose=FALSE, xlsxExport=TRUE,
                          fcTolerance=0.1){
  
  message("This is TPP version ", packageVersion("TPP"),".")
  
  ## ---------------------------------------------------------------------------
  ## 1) Import data and filter out rows for which column 'nonZeroCols'==0:
  datIn <- tppccrImport(configTable=configTable, data=data, idVar=idVar, 
                        fcStr=fcStr, naStrs=naStrs, qualColName=qualColName, 
                        nonZeroCols=nonZeroCols)
  expNames <- names(datIn)
  expNum   <- length(expNames)
  
  ## Extract directory from the filenames in config table, if specified:
  confgFields <- suppressMessages(importCheckConfigTable(infoTable=configTable, 
                                                         type="CCR"))
  files      <- confgFields$files
  outDirList <- importFct_makeOutputDirs(outDir=resultPath, fNames=files)
  flagDoWrite <- outDirList$doWrite
  pathDataObj <- outDirList$pathDataObj
  pathExcel=resultPath <- outDirList$outDir
  if (!flagDoWrite) plotCurves <- FALSE  
  
  ## ---------------------------------------------------------------------------
  ## 2) Normalize fold changes by their median:
  if (normalize){
    datNorm <- tppccrNormalize(data=datIn)
  } else {
    datNorm <- datIn
  }
  
  ## ---------------------------------------------------------------------------
  ## 3) Normalize to lowest concentrations and transform fold changes to [0,1]
  datNormToRef   <- tppccrNormalizeToReference(data=datNorm)
  datTransformed <- tppccrTransform(data=datNormToRef, fcCutoff=fcCutoff, 
                                    fcTolerance=fcTolerance)
  
  ## ---------------------------------------------------------------------------
  ## 4) calculate pEC50 values
  datFitted <- tppccrCurveFit(data=datTransformed, slopeBounds=slopeBounds, 
                              verbose=verbose, nCores=nCores)
  
  ## 5) Plot dose response curves
  if (plotCurves){
    datFitted <- tppccrPlotCurves(data=datFitted, resultPath=resultPath, 
                                  verbose=verbose, nCores=nCores,
                                  ggplotTheme=ggplotTheme)
  }
  
  ## 6) Produce results table
  resultTable <- tppccrResultTable(data=datFitted, r2Cutoff=r2Cutoff)
  
  ## 7) Save result table as data frame:
  if (flagDoWrite){
    save(list=c("resultTable"), file=file.path(pathDataObj, "results_TPP_CCR.RData"))
  }
  
  ## 8) Save result table as xlsx spreadsheet:
  if (xlsxExport){
    if (flagDoWrite){
      if (expNum > 1){
        expColors <- plotColors(expConditions=rep(NA, expNum), comparisonNums=rep(NA,expNum))
      } else {
        expColors <- NULL
      }
      tppExport(tab=resultTable, file=file.path(pathExcel, "results_TPP_CCR.xlsx"), expNames=expNames, expColors=expColors)
    } else {
      message("Cannot produce xlsx output because no result path is specified.")
    }
  }
  
  invisible(resultTable)
} 
