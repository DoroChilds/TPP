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
#' normData_longTable <- tpptrTidyUpESets(normData_eSets)$proteinMeasurements
#' hdacSubset <- subset(normData_longTable, grepl("HDAC", uniqueID))
#' hdacSplineFits <- tpptrFitSplines(data = hdacSubset, 
#'                                   factorsH1 = c("condition"))
#' 
#' @param data the data to be fitted
#' @param splineDF degrees of freedom for natural spline fitting.
#' @param factorsH1 which factors should be included in the alternative model?
#' @param factorsH0 which factors should be included in the null model?
#' @param computeAUC should areas under the spline curves be computed? 
#' Activation increases runtime requirements.
#' @param returnModels should the linear models be returned in a column of the
#' result table? Activation increases memory requirements.
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

tpptrFitSplines <- function(data, factorsH1, factorsH0 = c(), 
                            splineDF = 4, 
                            computeAUC = FALSE, returnModels = TRUE){
  ## Prepare data:
  dataGrouped <- data %>% group_by(uniqueID) 
  
  ## Define formulas for model fit:
  factorStrH1 <- paste0("factor(",factorsH1, ")", collapse = " * ") %>% paste0(" * ", .)
  factorStrH0 <- ifelse(length(factorsH0) > 0,
                        yes = paste0("factor(",factorsH0, ")", collapse = " * ") %>% 
                          paste0(" * ", .),
                        no = "")
  strFitH0 <- paste0("y ~ ns(x, df = ", splineDF, ")", factorStrH0)
  strFitH1 <- paste0("y ~ ns(x, df = ", splineDF, ")", factorStrH1)
  fitEqH0 <- as.formula(strFitH0)
  fitEqH1 <- as.formula(strFitH1)
  
  ## Fit null and alternative Model
  nProt <- length(unique(data$uniqueID))
  message("Fitting null models to ", nProt, " proteins.")
  fitsH0 <- dataGrouped %>% 
    do(fittedModel = fit_spline_model(dat = ., formula = fitEqH0))
  message("Fitting alternative models to ", nProt, " proteins.")
  fitsH1 <- dataGrouped %>% do(fittedModel = fit_spline_model(dat = ., formula = fitEqH1))
  
  # Combine results
  fitsCombined <- fitsH0 %>% mutate(testHypothesis = "null") %>%
    rbind(fitsH1 %>% mutate(testHypothesis = "alternative")) %>%
    mutate(splineDF = splineDF, testHypothesis = factor(testHypothesis)) %>% 
    ## Mark proteins where model fit was not successful
    ungroup %>%
    mutate(successfulFit = (sapply(fittedModel, class) != "try-error"))
  
  ## Evaluate goodness of fit of null and alternative models
  message("Evaluate goodness of fit of null and alternative models.")
  fitStats <- fitsCombined %>% 
    filter(successfulFit) %>%
    group_by(uniqueID, testHypothesis) %>%
    do(eval_spline_model(.$fittedModel[[1]]))
  
  if (computeAUC){
    aucTable <- tpptrSplineAUCs(data = data, 
                                 splineFits = fitsCombined, 
                                 factorsH1 = factorsH1)
    fitsCombined <- fitsCombined %>% 
      left_join(aucTable, by = c("uniqueID", "testHypothesis"))
  }
  
  if (!returnModels) fitsCombined <- fitsCombined %>% select(-fittedModel)
  
  ## Join fit stat df and model df (also include proteins for which model fit
  ## produced 'try-error')
  out <- fitsCombined %>% 
    left_join(fitStats, by = c("uniqueID", "testHypothesis"))
  
  return(out %>% arrange(uniqueID))
}

fit_spline_model <- function(dat, formula, algorithm = "lm"){
  if (algorithm == "lm"){
    fitResult <- try(lm(formula, data = dat), silent = TRUE)
  } else if (algorithm == "rlm"){ # Current algorithm for TPP-2D data fitting -> to do: re-use for TPP-2D fits.
    fitResult <- try(rlm(formula, data = dat, maxit = 50), silent = TRUE)
  }
  return(fitResult)
}

eval_spline_model <- function(lmFitObj){
  X <- model.matrix(lmFitObj)
  if (inherits(lmFitObj, "lmerMod")){
    nCoeffsRandom <- length(unlist(ranef(lmFitObj)))
    nCoeffsFixed <- length((fixef(lmFitObj)))
    nCoeffs <- nCoeffsRandom + nCoeffsFixed
  } else if (inherits(lmFitObj, "lm")){
    nCoeffs <- ncol(X)
  }
  fitStats <- data.frame(rss = sum(residuals(lmFitObj)^2), 
                         nCoeffs = nCoeffs, 
                         nObs = nrow(X), # protein numbers per group
                         sigma = sigma(lmFitObj), 
                         aicc = AICc(lmFitObj),
                         loglik = as.numeric(logLik(lmFitObj)))
  return(fitStats)
}
