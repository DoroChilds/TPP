#' @title Fit splines and perform f-Test
#'   
#' @description Fit splines through TR reference dataset and extrapolates relative 2D-TPP datapoints, 
#'  then compares spline fits of different treatments with non-treatment with an f-test 
#'  
#' @param data_2D result data.frame from a 2D-TPP CCR analysis
#' @param idVar character string indicating name of the columns containing the unique protein 
#'   identifiers
#' @param trRefDataPath character string with a link to a TPP-TR reference object RData file
#' @param fcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{fcStr} will be pasted in front of
#'   the names of the experiments.
#' @param refFcStr same as argument fcStr, but for the reference data.
#' @param resultPath location where to store dose-response curve plots and 
#'   results table.
#' @param ggplotTheme ggplot theme for protein-wise plots.
#' @param doPlot boolan value indicating whether protein-wise plots should be 
#'   produced Deactivating plotting decreases runtime.
#' @param verbose print description of problems for each protein for which splines fits could 
#'   not be performed
#'   
#' @examples 
#' data(panobinostat_2DTPP_smallExample)
#' config_tpp2d <- panobinostat_2DTPP_config
#' data_tpp2d <- panobinostat_2DTPP_data
#' trRef <- file.path(system.file("data", package="TPP"), "HepG2_trRefData.RData")
#' data2d <- tpp2dImportData(configTable = config_tpp2d,
#'  data = data_tpp2d, fcStr = NULL)
#' fcData2d <- tpp2dComputeFoldChanges(configTable = config_tpp2d,
#'                                     dataTable = data2d,
#'                                     intensityStr="sumionarea_protein_")
#' normData2d <- tpp2dDoMedianNorm(configTable = config_tpp2d,
#'                                 dataTable = fcData2d)
#' analysisResults <- tpp2dSplineFitAndTest(data_2D = normData2d,
#'                                          trRefDataPath = trRef,
#'                                          fcStr = "norm_rel_fc_protein_",
#'                                          refFcStr = "norm_rel_fc_protein_",
#'                                          doPlot = FALSE)
#'   
#' @return None
#' 
#' @export
tpp2dSplineFitAndTest <- function(data_2D, trRefDataPath, 
                                  idVar = "representative",
                                  fcStr = "norm_rel_fc_protein_", 
                                  refFcStr = "norm_rel_fc_protein_",
                                  resultPath = NULL,
                                  ggplotTheme = tppDefaultTheme(),
                                  doPlot = TRUE,
                                  verbose = FALSE){  
  
  # create tidy normalized 2D-TPP data
  dataList <- inferApparentStabilities(data_2D = data_2D, 
                                       trRefDataPath = trRefDataPath, 
                                       idVar = idVar,
                                       fcStr = fcStr,
                                       refFcStr = refFcStr)
  
  tppData_long_normalized <- dataList$normResult
  refTableLong <- dataList$referenceDataUsed
  
  if (!is.null(resultPath) && doPlot){
    if (!file.exists(file.path(resultPath))){
      dir.create(file.path(resultPath, plotDir), recursive=TRUE)
    }
    
    # plot splines
    splinePlots <- plotSplines(tppData_long_normalized, refTableLong)
    
    # export spline plots
    message("Exporting plots...")
    tpp2dExportPlots(plotList = splinePlots, outPath = resultPath, type = "spline")
    message("Done.")
  }else{
    message("Spline plots could not be saved, as no valid resultPath was indicated!")
  }
  
  # prepare for F test
  tppData_renamed <- dplyr::rename(tppData_long_normalized, 
                                   uniqueID = Protein_ID, 
                                   x = temperature, 
                                   y = fcNormalized)
  
  # do F test
  lmTable <- tpptrFitSplines(data=tppData_renamed, 
                             factorsH1=c("drugConc"), 
                             splineDF = 4)
  
  testResults <- tpptrFTest(fittedModels = lmTable, 
                            resultPath = NULL, 
                            doPlot = FALSE)
  
  # join results back to original table
  outTable <- left_join(as_data_frame(data_2D), testResults, by=setNames("uniqueID", idVar))
  
  return(as.data.frame(outTable))
}

