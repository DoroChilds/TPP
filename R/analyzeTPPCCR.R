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
#' @details Invokes the following steps: \enumerate{ \item Import data using the
#'   \code{\link{tppccrImport}} function. \item Perform normalization by fold 
#'   change medians (optional) using the \code{\link{tppccrNormalize}} function.
#'   To perform normalization, set argument \code{normalize=TRUE}. \item Fit and
#'   analyse dose response curves using the \code{\link{tppccrCurveFit}} 
#'   function. \item Export results to Excel using the \code{\link{tppExport}} 
#'   function. }
#'   
#'   The default settings are tailored towards the output of the python package 
#'   isobarQuant, but can be customised to your own dataset by the arguments 
#'   \code{idVar, fcStr, naStrs, qualColName}.
#'   
#'   If \code{resultPath} is not specified, the location of the input file 
#'   specified in \code{configTable} will be used. If the input data are not 
#'   specified in \code{configTable}, no result path will be set. This means 
#'   that no output files or dose response curve plots are produced and 
#'   \code{analyzeTPPCCR} just returns the results as a data frame.
#'   
#'   The function \code{analyzeTPPCCR} reports intermediate results to the 
#'   command line. To suppress this, use \code{\link{suppressMessages}}.
#'   
#' @details The dose response curve plots will be stored in a subfolder with 
#'   name \code{DoseResponse_Curves} at the location specified by 
#'   \code{resultPath}.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrResults <- analyzeTPPCCR(configTable=hdacCCR_config_repl1, 
#'                                data=hdacCCR_data_repl1)
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
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param normalize perform median normalization (default: TRUE).
#' @param ggplotTheme ggplot theme for dose response curve plots.
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param plotCurves boolan value indicating whether dose response curves should
#'   be plotted. Deactivating plotting decreases runtime.
#'   
#' @seealso tppDefaultTheme
#'   
#' @export
analyzeTPPCCR <- function(configTable, data=NULL, resultPath=NULL, 
                          idVar="gene_name", fcStr="rel_fc_", 
                          naStrs=c("NA", "n/d", "NaN", "<NA>"), qualColName="qupm",
                          normalize=TRUE, ggplotTheme=tppDefaultTheme(), 
                          nonZeroCols="qssm", 
                          r2Cutoff=0.8,  fcCutoff=1.5, slopeBounds=c(1,50),
                          plotCurves=TRUE){
  
  ## ---------------------------------------------------------------------------
  ## 1) Import data and filter out rows for which column 'nonZeroCols'==0:
  eSets <- tppccrImport(configTable=configTable, data=data, idVar=idVar, 
                        fcStr=fcStr, qualColName=qualColName, naStrs=naStrs, nonZeroCols=nonZeroCols)
  
  
  ## Extract directory from the filenames in config table, if specified:
  confgTableTmp <- importCheckConfigTable(infoTable=configTable)
  files        <- confgTableTmp$Path
  if (is.null(resultPath) & !is.null(files)){
    resultPath <- dirname(files[1])
    resultPath <- file.path(resultPath, "TPP_results")
  }
  if (is.null(resultPath)){
    message("No result directory specified. No output files or melting curve plots will be produced.") 
    plotCurves <- FALSE
  } else {
    message("Results will be written to ", resultPath)
    if (!file.exists(resultPath)) dir.create(resultPath, recursive=TRUE)
  }
  
  ## ---------------------------------------------------------------------------
  ## 2) Normalize fold changes by their median:
  if (normalize){
    eSetsNormalized <- tppccrNormalize(data=eSets)    
  } else {
    eSetsNormalized <- eSets
  }
  
  ## ---------------------------------------------------------------------------
  ## 3) Transform interesting proteins:
  ## 1. select 'stabilized' and 'destabilized' proteins
  ## 2. transform selected proteins
  eSetsFiltered <- tppccrTransform(data=eSetsNormalized, fcCutoff=fcCutoff)
    
  ## ---------------------------------------------------------------------------
  ## 4) calculate pEC50 values
  resultTable <- tppccrCurveFit(data=eSetsFiltered, resultPath=resultPath, 
                                ggplotTheme=ggplotTheme, doPlot=plotCurves, 
                                fcCutoff=fcCutoff, r2Cutoff=r2Cutoff, slopeBounds=slopeBounds)
  if (!is.null(resultPath)){
    ## 5) Save result table as data frame and Excel spreadsheet:
    save(list=c("resultTable"), file=file.path(resultPath, "results_TPP_CCR.RData"))
    tppExport(tab=resultTable, file=file.path(resultPath, "results_TPP_CCR.xlsx"))
  }
  
  
  return(resultTable)
} 
