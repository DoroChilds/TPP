#' @title Perform spline fitting
#' 
#' @description Fit natural splines to all proteins in a dataset.
#'   
#' @return A table containing the fitted models per protein
#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
#' normResults <- tpptrNormalize(data = tpptrData, 
#'                normReqs = tpptrDefaultNormReqs())
#' normData_eSets <- normResults$normData
#' normData_longTable <- tpptrTidyUpESets(normData_eSets)
#' hdacSubset <- subset(normData_longTable, grepl("HDAC", uniqueID))
#' hdacSplineFits <- tpptrFitSplines(data = hdacSubset, 
#'                                   factorsH1 = c("condition"), 
#'                                   nCores = 1)
#' 
#' @param data the data to be fitted
#' @param splineDF degrees of freedom for natural spline fitting.
#' @param factorsH1 which factors should be included in the alternative model?
#' @param factorsH0 which factors should be included in the null model?
#' @param computeAUC DEPRECATED
#' @param returnModels should the linear models be returned in a column of the
#' result table? Activation increases memory requirements.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' 
#' Argument \code{splineDF} specifies the degrees of freedom for natural spline 
#' fitting. As a single numeric value, it is directly passed on to the \code{splineDF} argument of 
#' \code{splines::ns}. Experience shows that \code{splineDF = 4} yields good results for TPP
#' data sets with 10 temperature points. It is also possible to provide a numeric vector. 
#' In this case, splines are fitted for each entry and the optimal value is chosen
#' per protein using Akaike's Information criterion.
#' 
#' @seealso \code{\link{ns}, \link{AICc}}
#' @export

tpptrFitSplines <- function(data, factorsH1, factorsH0 = character(0), 
                            splineDF = 3:7, 
                            computeAUC = NULL, returnModels = TRUE, 
                            nCores = "max"){
  
  ## ----------------------------------------------------------------------- ##
  ## General checks and preparation
  ## ----------------------------------------------------------------------- ##
  if (!missing(computeAUC)) warning("`computeAUC` is deprecated", call. = TRUE)
  
  if (!("uniqueID" %in% colnames(data)))
    stop("'data' must contain a column called 'uniqueID'")
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data", "factorsH1"))
  
  # ## Initialize variables to prevent "no visible binding for global
  # ## variable" NOTE by R CMD check:
  uniqueID = df = testHypothesis = fittedModel = successfulFit <- NULL
  
  # Produce informative error messages, if necessary:
  if (!is.character(factorsH0)) stop("'factorsH0' is of class ", class(factorsH0), ", but must be an atomic vector of class 'character'")
  if (!is.character(factorsH1)) stop("'factorsH1' is of class ", class(factorsH1), ", but must be an atomic vector of class 'character'")
  if (!is.numeric(splineDF)) stop("'splineDF' is of class ", class(splineDF), ", but must be an atomic vector of class 'numeric'")
  if (!all(factorsH1 %in% colnames(data))) stop("The column(s) '", paste(setdiff(factorsH1, colnames(data)), collapse = "', '"), "' specified by argument 'factorsH1' are not found in 'data'.")
  
  factorLevels <- subset(data, select = factorsH1) %>% purrr::map(. %>% unique)
  factorNums <- factorLevels  %>% purrr::map(. %>% length) %>% unlist
  factorsNonNA <- factorLevels  %>% purrr::map(. %>% is.na %>% all) %>% unlist
  if(any(factorsNonNA)) stop("At least one of the data column(s) specified by argument 'factorsH1' contain only NAs.") 
  if (any(factorNums <= 1)) stop("All of the data column(s) specified by arguments 'factorsH1' and 'factorsH0' need to contain varying levels.")
  
  ## ----------------------------------------------------------------------- ##
  ## Fit models for different degrees of freedom
  ## ----------------------------------------------------------------------- ##
  
  ## Define formulas for model fit:
  factorStrH1 <- paste0("factor(",factorsH1, ")", collapse = " * ") %>% paste0(" * ", .)
  factorStrH0 <- ifelse(length(factorsH0) > 0,
                        yes = paste0("factor(",factorsH0, ")", collapse = " * ") %>% paste0(" * ", .),
                        no = "")
  message(paste("Fitting smoothing splines and AICc values for the following degrees of freedom:", 
          paste(splineDF, collapse = ", ")))
  
  ## Loop over different degrees of freedom in parallel:
  nCores <- checkCPUs(cpus=nCores)
  doParallel::registerDoParallel(cores=nCores)
  t1 <- Sys.time()
  aicc_per_df <- foreach (df = splineDF, .combine=rbind) %dopar% {
    
    aiccCombined <- fit_splines_under_H0_and_H1(data = data, df = df,
                                                strH0 = factorStrH0, 
                                                strH1 = factorStrH1, 
                                                returnModels = FALSE)
    
    return(aiccCombined)
  }
  stopImplicitCluster() 
  
  if (!any(aicc_per_df$successfulFit)){
    stop("Spline smoothing did not converge for any protein. Consider using different degrees of freedom (parameter 'splineDF')")
  }
  
  ## ----------------------------------------------------------------------- ##
  ## Select desired model complexity individually per protein. 
  ## Criterion: corrected Aikake's information criterion (AICc) 
  ## ----------------------------------------------------------------------- ##
  message(paste("Select and re-fit models for selected degrees of freedom."))
  
  ## Select degrees of freedom that minimize the AICc:
  selectedDf <- modelSelector(fitStats = aicc_per_df %>% filter(successfulFit), 
                              criterion = "aicc",
                              hypothesis = "alternative") %>%
    right_join(aicc_per_df %>% distinct(uniqueID), by = "uniqueID")
  
  ## Re-fit models using the selected complexity:
  selectedModels <- selectedDf %>% filter(!is.na(splineDF)) %>%
    group_by(splineDF) %>%
    do({ 
      dfTmp <- unique(.$splineDF)
      idsTmp <- .$uniqueID %>% as.character()
      dataTmp <- data %>% filter(uniqueID %in% idsTmp)
      fit_splines_under_H0_and_H1(data = dataTmp, df = dfTmp,
                                  strH0 = factorStrH0,
                                  strH1 = factorStrH1,
                                  returnModels = TRUE)
    }) %>% right_join(aicc_per_df %>%distinct(uniqueID, testHypothesis, successfulFit), 
                      by = c("uniqueID", "testHypothesis", "successfulFit"))
  
  ## Evaluate goodness of fit of null and alternative models:
  message("Evaluate goodness of fit of null and alternative models.")
  fitStats <- selectedModels %>% 
    filter(successfulFit) %>%
    group_by(uniqueID, testHypothesis) %>%
    do(eval_spline_model(.$fittedModel[[1]])) %>% 
    arrange(uniqueID)
  
  ## ----------------------------------------------------------------------- ##
  ## Prepare output
  ## ----------------------------------------------------------------------- ##
  if (!returnModels) selectedModels <- selectedModels %>% select(-fittedModel)
  
  ## Join fit statistics and models (also include proteins for which model fit
  ## produced 'try-error')
  out <- selectedModels %>% 
    left_join(fitStats, by = c("uniqueID", "testHypothesis", "aicc")) %>%
    arrange(uniqueID)
   
  timeDiff <- Sys.time()-t1
  message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 2), " ", 
          units(timeDiff), "\n")
  
  return(out)
}
