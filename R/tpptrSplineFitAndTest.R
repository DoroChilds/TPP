#' @title Perform spline fitting and analyse by moderated F-test
#'   
#' @description A wrapper function around the functions \code{tpptrFitSplines}, 
#' \code{tpptrFTest}, \code{tpptrPlotSplines}, which fits natural splines to
#'   all proteins in a dataset and detect differential behaviour between
#'   conditions by a moderated F-test. The results are formatted as a wide table
#'   with one row per protein. This table contains all the original data, the
#'   test results, and (optionally) additional annotation columns for each 
#'   protein.
#'   
#' @return A list of two data frames: 1. A long table containing the spline
#'   predictions per protein and TMT-label 2. A long table containing the
#'   hypothesis test results per protein.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' normResults <- tpptrNormalize(data = tpptrData, 
#'                               normReqs = tpptrDefaultNormReqs())
#' normData_eSets <- normResults$normData
#' longTables <- normData_eSets %>% tpptrTidyUpESets
#' fitData <- longTables %>% extract2("proteinMeasurements")
#' proteinInfos <- longTables %>% extract2("proteinAnnotation")
#' hdacSplineFits <- tpptrSplineFitAndTest(data = fitData,
#'                                         factorsH1 = "condition",
#'                                         additionalCols = proteinInfos, 
#'                                         nCores = 1)
#' # Show estimated splines for HDAC1:
#' filter(hdacSplineFits, Protein_ID == "HDAC1")
#' # Quality control: test for replicate-specific effects:
#'  testResults <- tpptrSplineFitAndTest(data = fitData,
#'                                      factorsH1 = "replicate")
#' # -> Which proteins showed significant replicate effects?
#' testResults %>% filter(p_adj_NPARC <= 0.01) %>% select(Protein_ID, p_adj_NPARC)
#' 
#' @param data the data to be fitted.
#' @param resultPath location where to store the spline plots per protein.
#' @param ggplotTheme ggplot theme for melting curve plots.
#' @param doPlot boolan value indicating whether melting curves should be 
#'   plotted, or whether just the curve parameters should be returned.
#' @param splineDF degrees of freedom for natural spline fitting.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param verbose plot name of each fitted protein to the command lin as a means
#'   of progress report.
#' @param additionalCols additional annotation per protein to append to the 
#' result table.
#' @param factorsH1 which factors should be included in the alternative model?
#' @param factorsH0 which factors should be included in the null model?
#' @details Plots of the natural spline fits will be stored in a subfolder with 
#'   name \code{Spline_Fits} at the location specified by \code{resultPath}.
#'   
#'   Argument \code{splineDF} specifies the degrees of freedom for natural
#'   spline fitting. As a single numeric value, it is directly passed on to the
#'   \code{splineDF} argument of \code{splines::ns}. Experience shows that
#'   \code{splineDF = 4} yields good results for TPP data sets with 10
#'   temperature points. It is also possible to provide a numeric vector. In
#'   this case, splines are fitted for each entry and the optimal value is
#'   chosen per protein using Akaike's Information criterion.
#'   
#' @seealso \code{\link{ns}, \link{AICc}, \link{tppDefaultTheme}}
#' @export
tpptrSplineFitAndTest <- function(data,
                                  factorsH1,
                                  factorsH0 = c(),
                                  resultPath = NULL,
                                  ggplotTheme = tppDefaultTheme(),
                                  doPlot = TRUE,
                                  nCores = 'max', splineDF = 4, 
                                  additionalCols = NULL, # additional columns to display in the final output table
                                  verbose = FALSE){
  
  # ## Check whether plotting is possible (result path specified?)
  doPlot <- doPlot && !is.null(resultPath)
  
  lmTable <- tpptrFitSplines(data = data, 
                             factorsH1 = factorsH1,
                             factorsH0 = factorsH0,
                             splineDF = splineDF)
  testResults <- tpptrFTest(fittedModels = lmTable, 
                            resultPath = resultPath, doPlot = doPlot)
  
  if (doPlot){
    tpptrPlotSplines(data = data,
                     factorsH1 = factorsH1,
                     factorsH0 = factorsH0,
                     fittedModels = lmTable, 
                     testResults = testResults,
                     resultPath = resultPath, 
                     ggplotTheme = tppDefaultTheme())
  }
  
  dataList <- prepare_NPARC_results_for_export(measurements = data, 
                                               proteinInfos = additionalCols)
  testResultsRenamed <- testResults %>% rename(Protein_ID = uniqueID)
  resultTable <- mergeOutputTables_TR(dataList = dataList, 
                                      pValDF = testResultsRenamed, 
                                      qualCheckDF = NULL)
  
  return(resultTable)
}

prepare_NPARC_results_for_export <- function(measurements, proteinInfos = NULL){
  # Convert the data obtained by the 'tpptrSplineFitAnd Test' function into a 
  # format that can be used by the 'mergeOutputTables_TR' function in order
  # to create the final result table.
  
  emptyDF <- data.frame(Protein_ID = measurements$uniqueID) %>% distinct
  
  # Convert measurements from long to wide:
  fcDF <- measurements %>% 
    arrange(uniqueID, experiment, colName) %>%
    mutate(newColName = paste(colName, experiment, sep = "_")) %>%
    mutate(newColName = factor(newColName, levels = unique(newColName))) %>%
    select(uniqueID, y, newColName) %>% 
    spread(newColName, y) %>% 
    rename(Protein_ID = uniqueID)
  
  # Convert further protein annotation from long to wide:
  if (!is.null(proteinInfos)){
    rmCols <- c("plot", 
                meltCurveParamNames(returnParNames = TRUE, 
                                    returnPerformanceInfo = TRUE))
    
    infoTabfiltered <- proteinInfos %>% 
      arrange(uniqueID, experiment, variable) %>%
      filter(!variable %in% rmCols)
    
    if (nrow(infoTabfiltered) > 0){
      otherAnnotDF  <- infoTabfiltered %>%
        mutate(newColName = paste(variable, experiment, sep = "_")) %>%
        mutate(newColName = factor(newColName, levels = unique(newColName))) %>%
        select(uniqueID, value, newColName) %>% 
        spread(newColName, value) %>% 
        rename(Protein_ID = uniqueID)
    } else {
      otherAnnotDF <- emptyDF
    }
  } else {
    otherAnnotDF <- emptyDF
  }
  
  ## Store in a list that can be used by the 'mergeOutputTables_TR' function.
  dataList <- list(fcDF = fcDF,
                   curveParDF = emptyDF,
                   plotCol = emptyDF,
                   presenceDF = emptyDF,
                   modelInfoDF = emptyDF,
                   otherAnnotDF = otherAnnotDF)
  
  return(dataList)
}


