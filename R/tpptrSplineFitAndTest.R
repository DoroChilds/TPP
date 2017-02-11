#' @title Perform spline fitting and analyze by moderated F-test
#'   
#' @description A wrapper function around the functions \code{tpptrFitSplines}, 
#' \code{tpptrFTest}, \code{tpptrPlotSplines}, which fits natural splines to
#'   all proteins in a dataset and detect differential behavior between
#'   conditions by a moderated F-test. The results are formatted as a wide table
#'   with one row per protein. This table contains all the original data, the
#'   test results, and (optionally) additional annotation columns for each 
#'   protein.
#'   
#' @return A data frame in wide format with one row per protein. It contains 
#' the smoothing spline parameters and F-test results obtained by comparing 
#' the null and alternative models.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' fitData <- tpptrTidyUpESets(tpptrData)
#' hdacSplineFits <- tpptrSplineFitAndTest(data = fitData,
#'                                         factorsH1 = "condition",
#'                                         nCores = 1,
#'                                         doPlot = FALSE)
#' # Show estimated splines for HDAC1:
#' filter(hdacSplineFits, Protein_ID == "HDAC1")
#' # -> Which proteins showed significant condition effects?
#' hdacSplineFits %>% filter(p_adj_NPARC <= 0.01) %>% select(Protein_ID, p_adj_NPARC)

#' # Quality control: test for replicate-specific effects:
#'  testResults <- tpptrSplineFitAndTest(data = fitData,
#'                                      factorsH1 = "replicate",
#'                                      nCores = 1,
#'                                      doPlot = FALSE)
#' # -> Which proteins showed significant replicate effects?
#' testResults %>% filter(p_adj_NPARC <= 0.01) %>% select(Protein_ID, p_adj_NPARC)
#' 
#' @param data the data to be fitted.
#' @param resultPath location where to store the spline plots per protein.
#' @param ggplotTheme DEPRECATED.
#' @param doPlot boolean value indicating whether melting curves should be 
#'   plotted, or whether just the curve parameters should be returned.
#' @param splineDF degrees of freedom for natural spline fitting.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param additionalCols additional annotation per protein to append to the 
#' result table.
#' @param factorsH1 which factors should be included in the alternative model?
#' @param factorsH0 which factors should be included in the null model?
#' @param verbose DEPRECATED
#' 
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
#' @seealso \code{\link{ns}, \link{AICc}}
#' @export
tpptrSplineFitAndTest <- function(data,
                                  factorsH1,
                                  factorsH0 = character(),
                                  resultPath = NULL,
                                  doPlot = TRUE,
                                  nCores = 'max', splineDF = 3:7, 
                                  additionalCols = NULL, verbose = NULL,
                                  ggplotTheme = NULL){ # additional columns to display in the final output table
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data", "factorsH1"))
  
  if (!missing(verbose)) 
    warning("`verbose` is deprecated", call. = TRUE)
  
  if (!missing(ggplotTheme)) 
    warning("`ggplotTheme` is deprecated", call. = TRUE)
  
  if (!("uniqueID" %in% colnames(data)))
    stop("'data' must contain a column called 'uniqueID'")
  
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = path <- NULL
  
  ## Check whether plotting is possible (result path specified?)
  doPlot <- doPlot && !is.null(resultPath)
  
  ## Initialize data frames for result table creation:
  resultTableElements <- splitTidyMeasurementsForExport(
    measurements = data, 
    proteinInfos = additionalCols
  )
  
  lmTable <- tpptrFitSplines(data = data, 
                             factorsH1 = factorsH1,
                             factorsH0 = factorsH0,
                             splineDF = splineDF,
                             nCores = nCores)
  
  testResults <- tpptrFTest(
    fittedModels = lmTable, 
    resultPath = resultPath, doPlot = doPlot
  )  
  
  if (doPlot){
    paths <- tpptrPlotSplines(data = data,
                              fittedModels = lmTable, 
                              testResults = testResults,
                              resultPath = resultPath,
                              individual = TRUE,
                              overview = FALSE, 
                              returnPlots = FALSE,
                              control = list(nCores = 1, # Do not parallelize until this is more memory efficiently implemented.
                                             maxRank = 500, 
                                             highlightBelow = 0.05)
    )
    
    resultTableElements$plotCol <- paths$individual %>% 
      rename(Protein_ID = uniqueID, splinefit_plot = path)
  }
  
  resultTable <- mergeOutputTables_TR(
    dataList = resultTableElements, 
    pValDF = testResults %>% rename(Protein_ID = uniqueID), 
    qualCheckDF = NULL
  )
  
  return(resultTable)
}



