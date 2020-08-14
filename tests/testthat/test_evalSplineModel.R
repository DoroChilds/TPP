# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 2:4, nCores = 1)
)

modelH0 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "null") %>%
              extract2("fittedModel"))[[1]]

modelH1 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "alternative") %>%
              extract2("fittedModel"))[[1]]


test_that(desc = "allOk_H0", code = {
  
  mIn <- modelH0
  
  out <- TPP:::eval_spline_model(lmFitObj = mIn)
  
  check1 <- all(colnames(out) == c( "rss", "nCoeffs", "nObs", "sigma", "aicc", "loglik"))
  check2 <- out$rss == sum(residuals(mIn)^2)
  check3 <- out$loglik == logLik(mIn)
  
  expect_true(check1 & check2 & check3)
})


test_that(desc = "allOk_H1", code = {
  
  mIn <- modelH1
  
  out <- TPP:::eval_spline_model(lmFitObj = mIn)
  
  check1 <- all(colnames(out) == c( "rss", "nCoeffs", "nObs", "sigma", "aicc", "loglik"))
  check2 <- out$rss == sum(residuals(mIn)^2)
  check3 <- out$loglik == logLik(mIn)
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="modelMissing", code={
  
  expect_error(TPP:::eval_spline_model())
  
})


test_that(desc = "modelNULL", code = {
  
  mIn <- NULL
  
  out <- TPP:::eval_spline_model(lmFitObj = mIn)
  
  check1 <- all(colnames(out) == c( "rss", "nCoeffs", "nObs", "sigma", "aicc", "loglik"))
  check2 <- all(is.na(out[1,]))
  
  expect_true(check1 & check2)
  
})

test_that(desc="modelFitError", code={
  
  mIn <- try(lm(y ~ x, data = data.frame(x = NA, y = NA)), silent = TRUE)
  
  out <- TPP:::eval_spline_model(lmFitObj = mIn)
  
  check1 <- all(colnames(out) == c( "rss", "nCoeffs", "nObs", "sigma", "aicc", "loglik"))
  check2 <- all(is.na(out[1,]))
  
  expect_true(check1 & check2)
  
})

