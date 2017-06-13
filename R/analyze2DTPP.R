#' @title Analyze a 2D-TPP experiment
#'   
#' @description Performs the whole analysis workflow for 2D-TPP experiment by invoking routines 
#'   for data import, data processing, fold change computation, median normalization, 
#'   TPP-CCR curve fitting, plotting and production of the result table.
#'  
#' @return A data frame in which the model results (slopes and pEC50 values) are 
#'  stored row-wise for each protein and administered temperatures.
#'   
#' @references Becher, I., Werner, T., Doce, C., Zaal, E. A., Berkers, C. R., T"ogel, I., 
#' Salzer, E., Bantscheff, M., Savitski, M. M. (2016) 
#' Thermal profiling reveals phenylalanine hydroxylase as an off-target of panobinostat.
#' Nature Chemical Biology, 12(11), 908–910.
#' 
#' 
#' @details Invokes the following steps: \enumerate{ \item Import data using the
#'   \code{\link{tpp2dImport}} function. \item Remove zero sumionarea values. 
#'   \item Compute fold changes from raw data (sumionarea) 
#'   \item Perform normalization by fold 
#'   change medians (optional) using the \code{\link{tpp2dNormalize}} function.
#'   To perform normalization, set argument \code{normalize=TRUE}.}
#'  \code{paletteName} specifies the color palette to be used by the \code{\link{brewer.pal}} 
#' function from the \code{RColorBrewer} package to assign a separate color to 
#' each concentration.
#' 
#' 
#' @examples 
#' data(panobinostat_2DTPP_smallExample)
#' config_tpp2d <- panobinostat_2DTPP_config
#' data_tpp2d <- panobinostat_2DTPP_data
#' tpp2dResults <- analyze2DTPP(configTable = config_tpp2d, 
#'                              data = data_tpp2d,
#'                              methods=c("doseResponse"),
#'                              createReport="none",
#'                              nCores=1,
#'                              idVar = "representative",
#'                              addCol = "clustername",
#'                              intensityStr = "sumionarea_protein_",
#'                              nonZeroCols = "qusm")
#' 
#' @param configTable dataframe, or character object with the path to a file, 
#'   that specifies important details of the 2D-TPP experiment. See Section 
#'   \code{details} for instructions how to create this object.
#' @param data single dataframe, containing fold change measurements and 
#'   additional annotation columns to be imported. Can be used instead of 
#'   specifying the file path in the \code{configTable} argument.
#' @param resultPath location where to store dose-response curve plots and 
#'   results table.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the prefix \code{fcStr} 
#'   will be regarded as containing fold change values. Only relevant if 
#'   \code{compFC = FALSE}.
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the prefix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param methods vector of character strings that indicate which methods should be used 
#'   for the analysis (default: c("doseResponse"), alternative: c("splineFit") or 
#'   c("doseResponse", "splineFit")) 
#' @param addCol character vector indicating which additional columns to include 
#'   from the input data 
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param compFc boolean flag which indicates whether to perform fold change computation 
#'   regarding reference column from sumionareas (default: TRUE)
#' @param normalize perform median normalization (default: TRUE).
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param xlsxExport produce results table in xlsx format and store at the 
#' location specified by the \code{resultPath} argument.
#' @param plotAll boolean value indicating whether all dose response curves should
#' be generated. Deactivating plotting decreases runtime.
#' @param plotAllR2 boolean value indicating whether all dose response curves which
#' fulfill the demanded criteria (Rsquared, maximum plateau) should be generated. 
#' Deactivating plotting decreases runtime.
#' @param plotSingle boolean value indicating whether all dose response curves which
#' fulfill the demanded criteria (Rsquared, maximum plateau) should be generated. 
#' Deactivating plotting decreases runtime.
#' @param fractAbund boolean variable, if set to TRUE additional information concerning
#' sumionarea fractional abundance and dmso1 vs. dmso2 of adjacent temperatures is 
#' added to the output table
#' @param addInfo boolean variable, if set to TRUE additional information on counts of
#' stabilization and destabilization of each protein is added to the output table
#' @param trRef character string containing a valid system path to a previously generated TPP-TR
#' reference object
#' @param refFcStr character string indicating which columns in the reference data set contain 
#' the fold change values
#' @param createReport character string indicating whether a markdown report should be created
#'  and which format it have (default: "html_document", alternative: "pdf_document" or "none")
#' @param paletteName color palette (see details).
#' @param configFile DEPRECATED

#' 
#' @export
analyze2DTPP <- function(configTable, 
                         data = NULL, 
                         resultPath = NULL, 
                         idVar = "gene_name", 
                         fcStr = NULL, 
                         intensityStr = "signal_sum_",   
                         naStrs = c("NA", "n/d", "NaN", "<NA>"), 
                         methods = "doseResponse",
                         qualColName = "qupm", 
                         compFc = TRUE, 
                         normalize = TRUE, 
                         addCol = NULL,
                         nCores = 1, 
                         nonZeroCols = "qssm",
                         fcTolerance = 0.1, 
                         r2Cutoff = 0.8,  
                         fcCutoff = 1.5, 
                         slopeBounds = c(1,50), 
                         fractAbund = FALSE, 
                         xlsxExport = TRUE, 
                         plotAll = FALSE,
                         plotAllR2 = FALSE, 
                         plotSingle = FALSE, 
                         trRef = NULL, 
                         refFcStr="norm_rel_fc_", 
                         addInfo = FALSE, 
                         createReport = "none",
                         paletteName = "Spectral",
                         configFile) {
  
  message("This is TPP version ", packageVersion("TPP"),".")
  
  if (!missing(configFile)){
    warning("`configFile` is deprecated. Use 'configTable' instead.", call. = FALSE)
    configTable <- configFile
  }
  
  # # Check for missing function arguments
  # checkFunctionArgs(match.call(), c("configTable"))
  
  # import data
  datIn <- tpp2dImport(configTable=configTable, data=data, 
                       idVar=idVar, addCol=addCol, intensityStr=intensityStr, 
                       qualColName=qualColName, nonZeroCols = nonZeroCols,
                       fcStr=fcStr)
  
  ## Extract directory from the filenames in config table, if specified:
  confgFields <- suppressMessages(
    importCheckConfigTable(infoTable=configTable, type="2D")
    )
  files      <- suppressWarnings(confgFields$Path)
  outDirList <- importFct_makeOutputDirs(outDir=resultPath, fNames=files)
  flagDoWrite <- outDirList$doWrite
  resultPath <- outDirList$outDir
  if (!flagDoWrite) xlsxExport = plotAll = plotAllR2 = plotSingle <- FALSE
  
  # compute fold changes if requested
  if (compFc){
    fcStr <- "rel_fc_protein_"
    datIn <- tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStr) 
    
  }
  
  # do median normalization of fold changes 
  if (normalize){
    NormData2d <- tpp2dNormalize(data = datIn)
    
  }else{
    NormData2d <- datIn
  }
  # Make sure the TPP-CCR routine uses the correct columns, when there was 
  # normalization before:
  fcStrUpdated <- attr(NormData2d, "importSettings")$fcStrNorm
  
  # filter out row with no quality information
  if (length(which(is.na(NormData2d[[qualColName[1]]])))!=0){ # to do: shift this code into the curve fitting function to make it available for curve fitting outside this wrapper
    NormData2d <- NormData2d[-which(is.na(NormData2d[[qualColName[1]]])),]
  }
  
  # filter out proteins that have duplicated unique_ID to prevent crash 
  if (length(which(duplicated(NormData2d$unique_ID)))!=0){
    message("There are duplicated proteins in your experimental conditions! These are filtered out to run this analysis!
Please check your data quality and consider pre-filtering!")
    NormData2d <- NormData2d[-which(duplicated(NormData2d$unique_ID)),]
  }
  
  # calculate fractioncal abundance per curve
  if (fractAbund){
    NormData2d <- tpp2dCalcFractAbundance(data = NormData2d)
  }
  
  if ("doseResponse" %in% methods){
    
    # run TPP-CCR
    analysisResults <- tpp2dCurveFit(data = NormData2d, 
                                     nCores = nCores, 
                                     r2Cutoff = r2Cutoff, 
                                     fcCutoff = fcCutoff, 
                                     slopeBounds = slopeBounds,
                                     fcTolerance = fcTolerance)
    
    if (plotAll){
      # generate joint plots for all proteins detected
      plotList <- tpp2dCreateDRplots(data = analysisResults, type = "all", 
                                     verbose = TRUE, paletteName = paletteName)
      # write output file with plots
      tpp2dExportPlots(plotList = plotList, resultPath = resultPath, type = "all")
    }
    if (plotAllR2){
      # generate joint plots for all proteins detected with sufficient R2
      plotGoodList <- tpp2dCreateDRplots(data = analysisResults, type = "good", 
                                         verbose = TRUE, paletteName = paletteName)
      # write output file with plots
      tpp2dExportPlots(plotList = plotGoodList, resultPath = resultPath, type = "good")
    }
    if (plotSingle){
      # generate single plots for all protein in each condition fitted with sufficient R2
      plotSingleList <- tpp2dCreateDRplots(data = analysisResults, type = "single", verbose = TRUE)
      # write output file with plots
      tpp2dExportPlots(plotList = plotSingleList, resultPath = resultPath, type = "single")
    }
  }else{
    analysisResults <- NormData2d
  }
  
  # do spline fit over tpp-tr reference
  if (("splineFit" %in% methods) && !is.null(trRef)){
    # do f-test for splines fit
    analysisResults <- tpp2dSplineFitAndTest(data = analysisResults, 
                                             dataRef = trRef, 
                                             refIDVar = "Protein_ID",
                                             refFcStr = "norm_rel_fc_",
                                             resultPath = resultPath,
                                             doPlot = TRUE,
                                             verbose = TRUE,
                                             nCores = nCores)
    
  }else if("splineFit" %in% methods){
    message("The spline fit and corresponding f-Test could not be performed, as no TPP-TR reference dataset was specified!
Please check the file path you have specified for trRef!")
  }
  
  
  
  # add additional information e.g. how often protein was stabilized/destabilized if desired
  if (addInfo){
    analysisResults <- tpp2dAddAdditionalInfo(data = analysisResults, 
                                              idVar = idVar)
  }
  
  # add TR reference columns to result table
  if (!is.null(trRef) && file.exists(trRef)){
    analysisResults <-tpp2dMerge2dRef(data = analysisResults, 
                                      trRef = trRef, idVar = idVar)
  }
  
  # export results
  if (!is.null(resultPath) & xlsxExport){
    addPlotColumns <- any(c(plotAll, plotAllR2, plotSingle, !is.null(trRef)))
    
    tpp2dExport(tab = analysisResults, 
                outPath = resultPath, 
                addCol = addCol, addPlotColumns = addPlotColumns,
                trRef = trRef) 
  }
  
  # create markdown report
  if ((createReport!="none") && !is.null(resultPath)){
    tpp2dCreateReport(resultPath = resultPath, 
                      configFile = configTable, 
                      normalize = normalize,
                      configTable = configTable, 
                      data = analysisResults, 
                      idVar = "Protein_ID", 
                      fcStr = fcStr, 
                      fcStrUpdated = fcStrUpdated, 
                      documentType = createReport,
                      intensityStr = intensityStr, 
                      addCol = addCol, 
                      fcTolerance = fcTolerance, 
                      r2Cutoff = r2Cutoff, 
                      fcCutoff = fcCutoff, 
                      slopeBounds = slopeBounds, 
                      trRef = trRef)
  }
  
  # return output table
  return(analysisResults)
}
