#' @title Plot spline fits per protein
#' @description Plot spline fits per protein
#' 
#' @return None
#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' normResults <- tpptrNormalize(data = tpptrData, normReqs = tpptrDefaultNormReqs())
#' normData_eSets <- normResults$normData
#' normData_longTable <- tpptrTidyUpESets(normData_eSets)$proteinMeasurements
#' splineFits <- tpptrFitSplines(data = normData_longTable,
#'                               factorsH1 = "condition", returnModels = TRUE)
#' testResults <- tpptrFTest(fittedModels = splineFits, doPlot = FALSE)
#' tpptrPlotSplines(data = normData_longTable, fittedModels = splineFits,
#'                  factorsH1 = "condition", factorsH0 = c(),
#'                  testResults = testResults, resultPath = getwd(),
#'                  plotRanked = TRUE, maxRank = 20)
#' 
#' @param data the data to be plotted.
#' @param factorsH1 which factors were included in the alternative model?
#'        (necessary for correct prediction)
#' @param factorsH0 which factors were included in the null model? 
#'        (necessary for correct prediction)
#' @param fittedModels long table of fitted models. 
#'        Output of \code{\link{tpptrFitSplines}}.
#' @param testResults long table of p-values per protein. 
#'        Output of \code{\link{tpptrFTest}}.
#' @param resultPath location where to store the spline plots per protein.
#' @param ggplotTheme ggplot theme for melting curve plots.
#' @param plotAlphabetical Generate a summary pdf with 20 plots per page in 
#' alphabetical order?
#' @param plotRanked Generate a summary pdf with 20 plots per page, ordered by p-value?
#' @param maxRank if \code{plotRanked = TRUE}, how many of the top hits should 
#' be plotted (default: 500)?
#' 
#' @details Plots of the natural spline fits will be stored in a subfolder with 
#' name \code{Spline_Fits} at the location specified by \code{resultPath}.
#' 
#' @seealso \code{\link{ns}, \link{AICc}, \link{tppDefaultTheme}, 
#' \link{tpptrFitSplines}, \link{tpptrFTest}}
#' @export
tpptrPlotSplines <- function(data, factorsH1, factorsH0 = c(), 
                             fittedModels, testResults,
                             resultPath, ggplotTheme = tppDefaultTheme(),
                             plotAlphabetical = FALSE, 
                             plotRanked = TRUE, maxRank = 500){
  theme_set(ggplotTheme)
  a = 0.05
  highlTxt <- paste("adjusted p-value <=", a) 
  
  ## Define output directory and create it, if necessary:
  plotDir <- "Spline_Fits"
  if (!file.exists(file.path(resultPath, plotDir))){
    dir.create(file.path(resultPath, plotDir), recursive=TRUE)
  }
  
  ## Predict values across the whole range of the independent variable 
  ## (avoids re-fitting by geom_smooth):
  compFactorDF <- data %>% select_(factorsH1) %>% distinct %>% 
    mutate_(colorCol = factorsH1) %>% mutate(colorCol = factor(colorCol)) 
  xVec <- unique(data$x)
  xRange <- range(xVec)
  xDat <- compFactorDF %>% group_by_(factorsH1) %>%
    do(data.frame(x = seq(xRange[1], xRange[2], length.out = 50))) %>%
    ungroup
  message("Predict values for plotting based on the null and alternative models.")
  
  modelPred <- fittedModels %>%
    filter(successfulFit) %>% 
    group_by(uniqueID, testHypothesis) %>%
    do(data.frame(xDat, y = predict(.$fittedModel[[1]], newdata = xDat))) %>%
    ungroup %>%
    left_join(compFactorDF, by = factorsH1)
  
  modelPredH0 <- filter(modelPred, testHypothesis == "null") %>% 
    select(-testHypothesis)
  modelPredH1 <- filter(modelPred, testHypothesis == "alternative") %>% 
    select(-testHypothesis)
  
  ## Prepare data:
  ## 1. Test results:
  testResPlot <- testResults %>% 
    mutate(uniqueID = as.character(uniqueID)) %>%
    select(uniqueID, p_NPARC, p_adj_NPARC) %>%
    # Create text for p-value annotation of each plot:
    mutate(pStr = format(p_NPARC, scientific = TRUE, digits = 2),
           p.adjStr = format(p_adj_NPARC, scientific = TRUE, digits = 2),
           textStr = paste0("p.adj=", p.adjStr)) %>%
    # Assign page numbers for plotting ranked by p-value:
    arrange(p_adj_NPARC) %>%
    mutate(pageByP = assignBins(rev(1:nrow(.)), 20, collapseSmallest = FALSE)) %>%
    mutate(pageByP = plyr::mapvalues(pageByP, unique(pageByP), rev(unique(pageByP)))) %>%
    # Assign page numbers for plotting in alphabetical order:
    mutate(uniqueID = factor(as.character(uniqueID))) %>% # Re-assign levels, otherwise sorting by arrange() could not give the real alphabetical order
    arrange(uniqueID) %>%
    mutate(pageByID = assignBins(rev(1:nrow(.)), 20, collapseSmallest = FALSE)) %>%
    mutate(pageByID = mapvalues(pageByID, unique(pageByID), rev(unique(pageByID))))
  
  ## 2. Original data:
  data <- left_join(data, compFactorDF, by = factorsH1)
  
  ## Plot top N proteins (ranked by adjusted p-values):
  if (plotRanked){
    fName <- paste0("splineFit_top", maxRank, ".pdf")
    plotPath <- file.path(resultPath, plotDir, fName)
    testResTmp <- testResPlot %>% filter(rank(p_adj_NPARC) <= maxRank) %>% 
      mutate(pageNum = pageByP) %>% arrange(p_adj_NPARC)
    message('Plot the top N hits (ranked by adjusted p-values):')
    signfIDs <- testResTmp %>% filter(p_adj_NPARC <= a) %>% extract2("uniqueID")
    t1 <- Sys.time()
    plot_splines_to_file(data, modelPredH0, modelPredH1, 
                         testResPlot = testResTmp, plotPath = plotPath,
                         highlightIDs = signfIDs, highlightTxt = highlTxt)
    timeDiff <- Sys.time() - t1
    message("Runtime: ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  }
  
  ## Plot all proteins (alphabetical order):
  if (plotAlphabetical){
    fName <- "splineFit_alphabetical.pdf"
    plotPath <- file.path(resultPath, plotDir, fName)
    testResTmp <- testResPlot %>% mutate(pageNum = pageByID)
    message('Plot all proteins (alphabetical order):')
    signfIDs <- testResTmp %>% filter(p_adj_NPARC <= a) %>% extract2("uniqueID")
    t1 <- Sys.time()
    plot_splines_to_file(data, modelPredH0, modelPredH1, 
                         testResPlot = testResTmp, plotPath = plotPath,
                         highlightIDs = signfIDs, highlightTxt = highlTxt)
    timeDiff <- Sys.time() - t1
    message("Runtime: ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  }
  
  #   if (!is.na(dfs)) {
  #     plotTable <- mutate(plotTable, best_df=dfs)
  #   }
  return(NULL)
}

plot_splines_to_file <- function(data, modelPredH0, modelPredH1, testResPlot, 
                                 plotPath, highlightIDs, highlightTxt){
  
  pageVec <- sort(unique(testResPlot$pageNum))
  message("Plotting spline fits for ", nrow(testResPlot), " identifiers to file '", plotPath, "' .")
  message("Each page will display 20 individual plots (", length(pageVec)," pages in total).")
  
  # Sort items the same way as they are in 'testResPlot' to keep the same order for plotting 
  # (i.e. when proteins were sorted by p-value)
  sortIdx <- testResPlot %>% select(uniqueID) %>% mutate(idx = 1:nrow(.))
  data2 <- left_join(data, sortIdx, by = "uniqueID") %>% arrange(idx)
  
  # Plot to file
  pdf(file = plotPath, width = 11.811, height = 7.87, useDingbats = FALSE)
  for (pageTmp in pageVec){
    message("Plotting page ", pageTmp)
    testResTmp <- filter(testResPlot, pageNum == pageTmp)
    idsTmp <- testResTmp$uniqueID %>% as.character
    datTmp <- filter(data2, uniqueID %in% idsTmp)
    predH0Tmp <- filter(modelPredH0, uniqueID %in% idsTmp)
    predH1Tmp <- filter(modelPredH1, uniqueID %in% idsTmp)
    p <- create_spline_plots(datTmp = datTmp, predH0Tmp = predH0Tmp, 
                             predH1Tmp = predH1Tmp, testResTmp = testResTmp, 
                             highlightIDs = highlightIDs, 
                             highlightTxt = highlightTxt)
    print(p)
  }
  dev.off()
}

create_spline_plots <- function(datTmp, predH0Tmp, predH1Tmp, testResTmp, 
                                highlightIDs = c(), highlightTxt = ""){
  datTmp <- datTmp %>% 
    mutate(uniqueID = factor(uniqueID, levels = unique(uniqueID)))
  testResTmp <- testResTmp %>% 
    mutate(uniqueID = factor(uniqueID, levels = unique(uniqueID)))
  
  p <- ggplot(data = datTmp, aes(y = y, x = x)) + 
    ylim(-0.1, 1.5) +
    ylab("Fraction non-denatured") + 
    xlab("Temperature [\U00B0 C]") +
    scale_shape_discrete(name = "") +
    scale_color_discrete(name = "")
  # scale_fill_continuous(guide = guide_legend()) +
  # theme(legend.position="bottom")
  
  if (length(datTmp$uniqueID) >1) {
    p <- p + facet_wrap(~ uniqueID, nrow = 4, ncol = 5) +
      theme(strip.text = element_text(size = 7))
  } else{
    p <- p + ggtitle(unique(datTmp$uniqueID))
  }
  try(p <- p + geom_point(aes(shape = replicate, color = colorCol), size = 1, 
                          na.rm = TRUE))
  try(p <- p + geom_line(data = predH0Tmp, size = 0.7, na.rm = TRUE), 
      silent = TRUE)
  try(p <- p + geom_line(data = predH1Tmp, size = 0.7, na.rm = TRUE, 
                         aes(color = colorCol)), silent = TRUE)
  try(p <- p + geom_label(data = testResTmp, vjust = "top", hjust = "right", 
                          label.size = 0.5, inherit.aes = FALSE,
                          aes(label = textStr, x = Inf, y = Inf), alpha = 0.1), 
      silent = TRUE)
  
  
  backgroundTable <- datTmp %>% 
    select(uniqueID) %>% 
    distinct %>% 
    mutate(highlight = uniqueID %in% highlightIDs)
  p <- p + geom_rect(data = backgroundTable, aes(fill = factor(highlight)), 
                     xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                     alpha = 0.1, inherit.aes = FALSE) +
    scale_fill_manual(highlightTxt, values = c("FALSE" = "white", "TRUE" = "green"))
  
  return(p)
}
