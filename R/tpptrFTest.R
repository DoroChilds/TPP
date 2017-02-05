#' @title Analyze spline fits to detect differential behavior over time
#' 
#' @description Analyze fitted natural spline models and look for 
#' differential behaviour between conditions by a moderated F-test.
#'   
#' @return A long table containing the hypothesis test results per protein.
#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' normResults <- tpptrNormalize(data = tpptrData, normReqs = tpptrDefaultNormReqs())
#' normData_eSets <- normResults$normData
#' fitData <- tpptrTidyUpESets(normData_eSets)
#' fits <- tpptrFitSplines(data = fitData, factorsH1 = "condition", nCores = 1)
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
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  testHypothesis = rssH0 = rssH1 = nObsH1 = nCoeffsH1 = sigmaH1 = 
    residual_df_H1 = posterior_var_H1 = nCoeffsH0 = df2 = prior_df_H1 =
    F_moderated = df1 = F_scaled = df2_moderated = p_NPARC = uniqueID =
    F_statistic = p_adj_NPARC = staType = staValue <- NULL
  
  if (!("uniqueID" %in% colnames(fittedModels)))
    stop("'fittedModels' must contain a column called 'uniqueID'")
  
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
  
  ## to do: shift QC plots to an extra function 'tpptrFTestQCplots' and deprecate the doPlot and resultPath arguments
}
