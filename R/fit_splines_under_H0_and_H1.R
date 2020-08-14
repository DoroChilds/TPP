fit_splines_under_H0_and_H1 <- function(data, df, strH0, strH1, returnModels, intercept = TRUE){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = testHypothesis = fittedModel = successfulFit = splineDF <- NULL
  
  dataGrouped <- data %>% group_by(uniqueID) 
  nProt <- length(unique(data$uniqueID))
  
  #Fit Null Model
  strFitH0 <- paste0("y ~ ns(x, df = ", df, ")", strH0)
  fitEqH0 <- as.formula(strFitH0)
  
  message(paste0("Fitting null models to ", nProt,
                 " proteins (using ", df, " degrees of freedom)"))
  
  fitsH0 <- dataGrouped %>% 
    do(fittedModel = fit_spline_model(dat = ., formula = fitEqH0, 
                                      algorithm = "lm"))
  
  #Fit Alternative Model
  strFitH1 <- paste0("y ~ ns(x, df = ", df, ")", strH1)
  fitEqH1 <- as.formula(strFitH1)
  
  message(paste0("Fitting alternative models to ", nProt,
                 " proteins (using ", df, " degrees of freedom)"))
  
  fitsH1 <- dataGrouped %>% 
    do(fittedModel = fit_spline_model(dat = ., formula = fitEqH1, 
                                      algorithm = "lm"))
  
  # Combine results
  fitsCombined <- fitsH0 %>% mutate(testHypothesis = "null") %>%
    rbind(fitsH1 %>% mutate(testHypothesis = "alternative")) %>%
    ungroup() %>% 
    mutate(splineDF = df, testHypothesis = as.factor(testHypothesis)) %>% 
    ungroup %>%
    ## Mark proteins where model fit was not successful
    mutate(successfulFit = (sapply(fittedModel, class) != "try-error"))
  
  aiccValues <- fitsCombined %>%
    filter(successfulFit) %>%
    group_by(uniqueID, splineDF, testHypothesis) %>%
    do(data.frame(aicc = AICc(.$fittedModel[[1]])))
  
  fitsCombined <- left_join(fitsCombined, aiccValues, 
                   by = c("uniqueID", "splineDF", "testHypothesis"))
  
  if (returnModels == FALSE){
    fitsCombined <- fitsCombined %>% select(-fittedModel)
  }
  
  return(fitsCombined)
}