#' @title Analyse spline fits to detect differential behaviour over time
#' 
#' @description Analyse fitted natural spline models and look for 
#' differential behaviour between conditions by a moderated F-test.
#'   
#' @return A long table containing the hypothesis test results per protein.
#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' normResults <- tpptrNormalize(data = tpptrData, normReqs = tpptrDefaultNormReqs())
#' normData_eSets <- normResults$normData
#' longTables <- normData_eSets %>% tpptrTidyUpESets
#' fitData <- longTables %>% extract2("proteinMeasurements")
#' fits <- tpptrFitSplines(data = fitData, factorsH1 = "condition")
#' testResults <- tpptrFTest(fittedModels = fits)
#' 
#' @param fittedModels a table of fitted spline models (produced by \code{tpptrFitSplines}).
#' @param doPlot boolean value indicating whether QC plots should be produced.
#' Currently, QC plots comprise distributions of the F statistics, and the 
#' p-values before/ after Benjamini Hochberg adjustment.
#' @param resultPath location where to store QC plots, if \code{doPlot} = TRUE.
#' 
#' @details If \code{doPlot} is \code{TRUE}, but no \code{resultPath} is 
#' specified, the plots will be prompted to the active device.
#' 
#' The moderated F-statistic is calculated by the following equation:
#' ...
#' 
#' @seealso \code{\link{ns}, \link{squeezeVar}}
#' @export
tpptrFTest <- function(fittedModels, doPlot = FALSE, resultPath = NULL){
  # results <- modelSelector(testResults=results, method=modMet, criterion=modCrit)
  
  ## Create wide table of goodness of fit statistics
  ## (not the most elegant solution, but spread only works for single column)
  fitStatsH0 <- fittedModels %>% filter(testHypothesis == "null") %>% 
    mutate(fittedModel = NULL) %>% # remove column "fittedModel", if it exists at all
    select( -testHypothesis)
  fitStatsH1 <- fittedModels %>% filter(testHypothesis == "alternative") %>% 
    mutate(fittedModel = NULL) %>% # remove column "fittedModel", if it exists at all
    select( -testHypothesis)
  fitStatsBothModels <- full_join(fitStatsH0, fitStatsH1, by = "uniqueID", 
                                  suffix = c("H0", "H1")) 
  
  testResults <- fitStatsBothModels %>% 
    ungroup %>%
    # Compute ordinary F-statistic as suggested by Storey et al. 2005:
    mutate(F_statistic = (rssH0-rssH1)/rssH1) %>% 
    # Compute moderated F-statistic:
    mutate(residual_df_H1 = nObsH1-nCoeffsH1,   
           posterior_var_H1 = squeezeVar(sigmaH1^2, residual_df_H1)$var.post, 
           prior_df_H1 = squeezeVar(sigmaH1^2, residual_df_H1)$df.prior,
           F_moderated = (rssH0-rssH1)/posterior_var_H1) %>%
    # Compute p-vals using F-distribution (requires scaling of mod. F-statistic, 
    # see Smyth (2004) p. 14 for details):
    mutate(df1 = nCoeffsH1-nCoeffsH0,         # numerator degrees of freedom
           df2 = (nObsH1-nCoeffsH1), # denominator degrees of freedom
           df2_moderated = df2 + prior_df_H1, # adjustment by adding prior estimator
           F_scaled = F_moderated*(1/df1),       # scale moderated F-statistic
           p_NPARC = 1-suppressWarnings(pf(F_scaled, df1, df2_moderated)),  # calculate p-val. Suppress warnings due to NAs
           p_adj_NPARC = p.adjust(p_NPARC, method = "fdr")) %>% # Benjamini-Hochberg correction
    # Select relevant columns (remove values already contained in input table)
    select(uniqueID, F_statistic, F_moderated, F_scaled, residual_df_H1, 
           prior_df_H1, df1, df2, df2_moderated, posterior_var_H1, 
           p_NPARC, p_adj_NPARC)
  
  if (doPlot){
    if(!is.null(resultPath)){
      pdf(file = file.path(resultPath, "QCplots_fTest.pdf"), width=8, height=9)
    }
    message("Creating QC plots to visualize test statistics and p-values")

    # Plot distributions of F-statistics
    plotList <- testResults %>% na.omit %>%
      gather(staType, staValue, c(F_statistic, F_moderated, F_scaled)) %>% 
      mutate(staType = factor(staType)) %>%
      group_by(staType) %>%
      do(plotObj = plot_fSta_distribution(dataLong = .))
    print(plotList$plotObj)

    # Plot distributions of p-values and adjusted p-values
    plotList <- testResults %>% na.omit %>%
      group_by(df1, df2, df2_moderated) %>%
      do(plotObj = plot_pVal_distribution(dataWide = .))
    print(plotList$plotObj)

    if(!is.null(resultPath)) dev.off()
    message("done.\n")
  }
  
  return(testResults)
}

plot_pVal_distribution <- function(dataWide){
  plotDat <- dataWide %>% 
    select(uniqueID, p_NPARC, p_adj_NPARC, df1, df2, df2_moderated) %>%
    gather(pValType, pValue, c(p_NPARC, p_adj_NPARC)) %>%
    mutate(pValType = factor(pValType), 
           df2_moderated = round(df2_moderated, 2)) %>%
    na.omit()
  
  mainTitle = "P-values before and after Benjamini-Hochberg adjustment"
  subTitle = paste("df1 = ", unique(plotDat$df1),
                 ", df2 = ", unique(plotDat$df2),
                 ", df2_moderated = ", unique(plotDat$df2_moderated), sep = "")
  numProt = plotDat %>% group_by(pValType)%>% summarize(n = n()) %>%
    mutate(label = paste("n =", n))
  
  p <- ggplot(data = plotDat, aes(x = pValue)) +
    geom_histogram(aes(y=..density../max(..density..)), # Histogram with density instead of count on y-axis
                   alpha = 1, binwidth = 0.05, na.rm = TRUE) +
    geom_text(data = numProt, aes(x = 0.1, y = 1, label = label), 
              color = "red", inherit.aes = FALSE) +
    ylab("Relative frequency") +
    ylim(c(-0.1,1.1)) + xlim(c(-0.1,1.1)) +
    # geom_density(alpha = .2, fill = "#FF6666") + # Overlay with transparent density plot
    facet_grid(pValType ~ .) +
    ggtitle(paste(mainTitle, subTitle, sep = "\n"))
  return(p)
}

plot_fSta_distribution <- function(dataLong){
  plotDat <- dataLong %>% 
    select(uniqueID, staType, staValue, df1, df2, df2_moderated) %>%
    mutate(staType = factor(staType), df2_moderated = round(df2_moderated, 2)) %>%
    na.omit()
  
  mainTitle = "F-statistics and moderated F-statistics"
  subTitle = paste("Type = ", unique(plotDat$staType), sep = "")
  numProt = plotDat %>% group_by(staType, df1, df2, df2_moderated)%>% 
    summarize(n = n(), DF1 = unique(df1), DF2 = unique(df2), MOD.DF2 = unique(df2_moderated)) %>%
    mutate(label = paste("n =", n, ", df1 =", DF1, ", df2 =", df2, ", df2_moderated =", df2_moderated, sep = " "))
  
  p <- ggplot(data = plotDat, aes(x = staValue)) +
    geom_histogram(aes(y=..density../max(..density..)), # Histogram with density instead of count on y-axis
                   alpha = 1, bins = min(nrow(plotDat), 100), na.rm = TRUE) +
    geom_text(data = numProt, aes(x = 0.1, y = 1, label = label), 
              color = "red", inherit.aes = FALSE, hjust = "left", vjust = "top") +
    # geom_density(alpha = .2, fill = "#FF6666") + # Overlay with transparent density plot
    facet_grid(df1 + df2 + df2_moderated ~ .) +
    ylim(c(0,1)) +
    ylab("Relative frequency") + xlab("F statistic") +
    ggtitle(paste(mainTitle, subTitle, sep = "\n"))
  # stat_function(fun = df, colour = "red", args = list(df1 = df1, df2 = df2))
  
  return(p)
}