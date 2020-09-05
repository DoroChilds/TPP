#' @title Fit splines and perform f-Test
#'   
#' @description Fit splines through TR reference dataset and extrapolates relative 2D-TPP datapoints, 
#'  then compares spline fits of different treatments with non-treatment with an f-test 
#' 
#' @param data_2D DEPRECATED
#' @param data result data.frame from a 2D-TPP CCR analysis
#' @param trRefDataPath DEPRECATED
#' @param dataRef reference data from a TPP TR analysis on the same cell line as
#' @param refIDVar character string indicating name of the columns containing the unique protein 
#'   identifiers in the reference data set
#' @param refFcStr character string indicating which columns contain the actual 
#'   fold change values in the reference data. The suffix \code{fcStr} will be 
#'   pasted in front of
#'   the names of the experiments.
#' @param resultPath location where to store dose-response curve plots and 
#'   results table.
#' @param ggplotTheme DEPRECATED
#' @param doPlot boolean value indicating whether protein-wise plots should be 
#'   produced Deactivating plotting decreases runtime.
#' @param verbose print description of problems for each protein for which splines fits could 
#'   not be performed
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#'   
#' @details 
#' dataRef can either be a tidy data frame of TPP-TR reference data,
#' a list with TPP-TR reference data and additional information produced by 
#' \code{\link{tpp2dCreateTPPTRreference}}, or a character string with a link to 
#' the data in one of the described formats.

#'   
#' @examples 
#' data(panobinostat_2DTPP_smallExample)
#' config_tpp2d <- panobinostat_2DTPP_config
#' data_tpp2d <- panobinostat_2DTPP_data
#' trRef <- file.path(system.file("data", package="TPP"), 
#'   "TPPTR_reference_results_HepG2.RData")
#' datIn <- tpp2dImport(configTable = config_tpp2d,
#'                       data = data_tpp2d,
#'                       idVar = "representative",
#'                       addCol = "clustername",
#'                       intensityStr = "sumionarea_protein_",
#'                       nonZeroCols = "qusm")
#' fcData2d <- tpp2dComputeFoldChanges(data = datIn)
#' normData2d <- tpp2dNormalize(data = fcData2d)
#' analysisResults <- tpp2dSplineFitAndTest(data = normData2d,
#'                                          dataRef = trRef,
#'                                          refIDVar = "Protein_ID",
#'                                          refFcStr = "norm_rel_fc_protein_",
#'                                          doPlot = FALSE,
#'                                          nCores = 1)
#'   
#' @return None
#' 
#' @export
tpp2dSplineFitAndTest <- function(data_2D = NULL, 
                                  data, 
                                  trRefDataPath = NULL, 
                                  dataRef, 
                                  refIDVar = "Protein_ID",
                                  refFcStr = "norm_rel_fc_",
                                  resultPath = NULL,
                                  doPlot = TRUE,
                                  verbose = FALSE,
                                  nCores = "max",
                                  ggplotTheme = NULL){ 
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  plotDir = temperature = fcNormalized = uniqueID <- NULL
  
  if (!missing(data_2D)){
    warning("`data_2D` is deprecated.", call. = TRUE)
  }
  
  if (!missing(trRefDataPath)){
    warning("`trRefDataPath` is deprecated.", call. = TRUE)
  }
  
  if (!missing(ggplotTheme)){
    warning("`ggplotTheme` is deprecated.", call. = TRUE)
  }
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data", "dataRef"))
  
  # create tidy normalized 2D-TPP data
  dataList <- inferApparentStabilities(data_2D = data, 
                                       dataRef = dataRef, 
                                       refIDVar = refIDVar,
                                       refFcStr = refFcStr)
  
  tppData_long_normalized <- dataList$normResult
  refTableLong <- dataList$referenceDataUsed
  
  if (doPlot){
    if (!is.null(resultPath)){
      if (!file.exists(file.path(resultPath))){
        dir.create(file.path(resultPath, plotDir), recursive=TRUE)
      }
      
      # plot splines
      splinePlots <- plot_2D_data_on_temperature_range(
        tppData_long_normalized, 
        refTableLong
        )
      
      # export spline plots
      message("Exporting plots...")
      tpp2dExportPlots(plotList = splinePlots, resultPath = resultPath, 
                       type = "spline")
      message("Done.")
    }else{
      message("Spline plots could not be saved, as no valid resultPath was 
              indicated!")
    }
  }
  
  
  # prepare for F test
  tppData_renamed <- dplyr::rename(tppData_long_normalized, 
                                   x = temperature, 
                                   y = fcNormalized)
  
  # do F test
  lmTable <- tpptrFitSplines(data=tppData_renamed, 
                             factorsH1 = c("drugConc"),
                             factorsH0 = character(0),
                             splineDF = 4,
                             returnModels = TRUE, 
                             nCores = nCores)
  

  testResults <- tpptrFTest(fittedModels = lmTable, 
                            resultPath = NULL, 
                            doPlot = FALSE)
  
  # join results back to original table
  idVar <- checkAndReturnDataSetting(attr(data, "importSettings"), 
                                     "proteinIdCol", colnames(data))
  
  outTable <- left_join(as_tibble(data), testResults, by=setNames("uniqueID", idVar))
  
  return(as.data.frame(outTable))
}

