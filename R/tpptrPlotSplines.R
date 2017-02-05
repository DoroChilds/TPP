#' @title Plot spline fits per protein
#' @description Plot spline fits per protein
#' 
#' @return None
#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' tidyData <- tpptrTidyUpESets(tpptrData)
#' splineFits <- tpptrFitSplines(data = tidyData, nCores = 1,
#'                               factorsH1 = "condition", returnModels = TRUE)
#' testResults <- tpptrFTest(fittedModels = splineFits, doPlot = FALSE)
#' tpptrPlotSplines(data = tidyData, fittedModels = splineFits,
#'                  individual = FALSE,
#'                  testResults = testResults, resultPath = getwd())
#' 
#' @param data long table of proteins measurements that were used for spline fitting.
#' @param factorsH1 DEPRECATED
#' @param factorsH0 DEPRECATED
#' @param fittedModels long table of fitted models. 
#'        Output of \code{\link{tpptrFitSplines}}.
#' @param testResults long table of p-values per protein. 
#'        Output of \code{\link{tpptrFTest}}.
#' @param resultPath an optional character vector with the name of the path where the plots should be saved.
#' @param individual logical. Export each plot to individual files?
#' @param overview logical. Generate summary pdfs?
#' @param returnPlots logical. Should the ggplot objects be returned as well?
#' @param control a list of general settings. 
#' @param maxRank DEPRECATED
#' @param highlightBelow DEPRECATED
#' @param plotIndividual DEPRECATED
#' @param plotAlphabetical DEPRECATED
#' Contains the following fields:
#' \itemize{ 
#' \item{\code{nCores}: number of CPUs for parallel production of plots per 
#' protein if \code{individual = TRUE} (default: "max")}
#' \item{\code{maxRank}: how many of the top hits should 
#' be plotted if \code{overview = TRUE} (default: 500)} 
#'  \item{\code{highlightBelow}: maximum adjusted p-value 
#' for which a protein is highlighted by a different background color if 
#' \code{overview = TRUE} (default: 0.05)}
#' }
#' 
#' @details Plots of the natural spline fits will be stored in a subfolder with 
#' name \code{Spline_Fits} at the location specified by \code{resultPath}.
#' 
#' Exporting each plot to individual files (individual = TRUE) can 
#' cost runtime and the resulting files can be tedious to browse. 
#' If you just want to browse the results, use \code{overview = TRUE} 
#' instead.
#' 
#' If \code{overview = TRUE}, two summary PDFs are created that enable quick 
#' browsing through all results.  They contain the plots in alphacetical order 
#' (\code{splineFit_alphabetical.pdf}), or ranked by p-values 
#' (\code{splineFit_top_xx.pdf}, where xx is the maximum rank defined by 
#' \code{overviewSettings$maxRank}). 
#' 
#' 
#' @seealso \code{\link{ns}, \link{AICc}, 
#' \link{tpptrFitSplines}, \link{tpptrFTest}}
#' @export
tpptrPlotSplines <- function(data, factorsH1 = NULL, factorsH0 = NULL, 
                             fittedModels, testResults,
                             resultPath = NULL,
                             individual = TRUE,
                             overview = FALSE,
                             returnPlots = FALSE,
                             control = list(nCores = "max",
                                            maxRank = 500, 
                                            highlightBelow = 0.05),
                             maxRank = NULL, highlightBelow = NULL,
                             plotIndividual = NULL, plotAlphabetical = NULL){
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), 
                    c("data", "fittedModels", "testResults"))
  
  if (!missing(factorsH1)) 
    warning("`factorsH1` is deprecated", call. = TRUE)
  
  if (!missing(factorsH0)) 
    warning("`factorsH0` is deprecated", call. = TRUE)
  
  if (!missing(maxRank)) 
    warning("`maxRank` is deprecated", call. = TRUE)
  
  if (!missing(highlightBelow)) 
    warning("`highlightBelow` is deprecated", call. = TRUE)
  
  if (!missing(plotIndividual)) 
    warning("`plotIndividual` is deprecated", call. = TRUE)
  
  if (!missing(plotAlphabetical)) 
    warning("`plotAlphabetical` is deprecated", call. = TRUE)
  
  if (!("uniqueID" %in% colnames(data)))
    stop("'data' must contain a column called 'uniqueID'")
  
  if (!("uniqueID" %in% colnames(fittedModels)))
    stop("'fittedModels' must contain a column called 'uniqueID'")
  
  if (!("uniqueID" %in% colnames(testResults)))
    stop("'testResults' must contain a column called 'uniqueID'")
  
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = p_NPARC = p_adj_NPARC = splineDF = p.adjStr = pageByP = pageByID =
    textString = path <- NULL
  
  ## If export to files is desired, define paths first
  doPlot <- !is.null(resultPath)
  
  if (doPlot){
    
    plotDir <- file.path(resultPath,"Spline_Fits")
    
  } else {
    
    plotDir <- NULL
    
    if (!returnPlots){
      warning("'resultPath' is set to NULL and 'returnPlots' is set to FALSE. With these settings, no output will be produced and the plots will be lost.", call. = TRUE)
    }
  }
  
  ## Prepare test results:
  annotatedTestResults <- testResults %>%
    select(uniqueID, p_NPARC, p_adj_NPARC) %>%
    left_join(distinct(fittedModels, uniqueID, splineDF), by = "uniqueID") %>%
    # Create text for p-value annotation of each plot:
    mutate(pStr = format(p_NPARC, scientific = TRUE, digits = 2),
           p.adjStr = format(p_adj_NPARC, scientific = TRUE, digits = 2),
           textString = paste0("p.adj=", p.adjStr,"\n", 
                               "spline degrees of freedom=",splineDF)) %>%
    # add empty rows for proteins with unsuccessful model fits
    right_join(distinct(data, uniqueID), by = "uniqueID") %>% 
    # convert factors to characters to enable sorting by p-values:
    mutate(uniqueID = as.character(uniqueID)) %>%
    # Assign page numbers for plotting ranked by p-value:
    arrange(p_adj_NPARC) %>%
    mutate(pageByP = assignBins(rev(1:nrow(.)), 20, collapseSmallest = FALSE)) %>%
    mutate(pageByP = plyr::mapvalues(pageByP, unique(pageByP), rev(unique(pageByP)))) %>%
    # Re-assign levels, otherwise sorting by arrange() could not give the real 
    # alphabetical order:
    mutate(uniqueID = factor(as.character(uniqueID))) %>% 
    # Assign page numbers for plotting in alphabetical order:
    arrange(uniqueID) %>%
    mutate(pageByID = assignBins(rev(1:nrow(.)), 20, collapseSmallest = FALSE)) %>%
    mutate(pageByID = mapvalues(pageByID, unique(pageByID), rev(unique(pageByID)))) %>%
    # Convert back to character:
    mutate(uniqueID = as.character(uniqueID))
  
  if(individual){
    ## Generate plots per protein:
    plotAnnotation = annotatedTestResults %>% select(uniqueID, textString)
    
    individualPlots <- plotIndividual(data = data, 
                                      fittedModels = fittedModels, 
                                      plotAnnotation = plotAnnotation, 
                                      plotDir = plotDir, 
                                      filePrefix = "smoothingSpline", 
                                      returnPlots = returnPlots,
                                      nCores = control$nCores)
    
    if (!is.null(resultPath)){
      individualPlots <- individualPlots %>%
        mutate(path = gsub(resultPath, "", path))
    }
    
  } else individualPlots <- NULL
  
  overviewPlots <- NULL
  
  # highlTxt <- paste("adjusted p-value <=", highlightBelow) 
  # 
  # ## Plot top N proteins (ranked by adjusted p-values):
  # if (overview){
  #   
  #   message('Plot the top ', maxRank, ' hits (ranked by adjusted p-values):')
  #   
  #   fName <- paste0("splineFit_top", maxRank, ".pdf")
  #   
  #   plotPath <- file.path(resultPath, plotDir, fName)
  #   
  #   testResTmp <- plotText %>% filter(rank(p_adj_NPARC) <= maxRank) %>% 
  #     mutate(pageNum = pageByP) %>% arrange(p_adj_NPARC)
  #   
  #   plotsTmp <- testResTmp %>%
  #     distinct(uniqueID, pageNum) %>%
  #     left_join(allPlots, by = "uniqueID") 
  #   
  #   signfIDs <- testResTmp %>% 
  #     filter(p_adj_NPARC <= highlightBelow) %>% 
  #     extract2("uniqueID")
  #   
  #   t1 <- Sys.time()
  #   
  #   plot_splines_to_file(plotsTmp, plotPath = plotPath,
  #                        highlightIDs = signfIDs, highlightTxt = highlTxt)
  #   timeDiff <- Sys.time() - t1
  #   message("Runtime: ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  # }
  # 
  # ## Plot all proteins (alphabetical order):
  # if (overview){
  #   fName <- "splineFit_alphabetical.pdf"
  #   plotPath <- file.path(resultPath, plotDir, fName)
  #   testResTmp <- plotText %>% mutate(pageNum = pageByID)
  #   message('Plot all proteins (alphabetical order):')
  #   signfIDs <- testResTmp %>% filter(p_adj_NPARC <= a) %>% extract2("uniqueID")
  #   t1 <- Sys.time()
  #   plot_splines_to_file(data, modelPredH0, modelPredH1, 
  #                        plotText = testResTmp, plotPath = plotPath,
  #                        highlightIDs = signfIDs, highlightTxt = highlTxt)
  #   timeDiff <- Sys.time() - t1
  #   message("Runtime: ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  # }
  # 
  #   if (!is.na(dfs)) {
  #     plotTable <- mutate(plotTable, best_df=dfs)
  #   }
  
  out <- list(individual = individualPlots, overview = overviewPlots)
  
  return(out)
}
