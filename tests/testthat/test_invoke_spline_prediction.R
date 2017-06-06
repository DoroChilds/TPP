# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))

splineFits <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 4, nCores = 1)
)

splineFits2 <- suppressMessages(
  tpptrFitSplines(data = testData,
                  factorsH1 = c("condition", "replicate"), 
                  factorsH0 = c("replicate"),
                  returnModels = TRUE, 
                  splineDF = 4,
                  nCores = 1)
)

xNew <- seq(0, 100, by = 10)

test_that(desc="allOk", code={
  
  fitsIn <- splineFits
  xIn <- xNew
  
  predictions <- TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn)
  
  check1 <- nrow(predictions) == (length(xIn) * (nrow(fitsIn) + sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit,na.rm = TRUE)))
  check2 <- setdiff(colnames(fitsIn), colnames(predictions)) == "fittedModel"
  check3 <- all(setdiff(colnames(predictions), colnames(fitsIn)) == c("x", "y", "condition"))

  expect_true(check1 & check2 & check3)
  
})

test_that(desc="allOk2", code={
  
  fitsIn <- splineFits2 %>% filter(uniqueID %in% c("HDAC1", "HDAC2"))
  xIn <- xNew
  
  predictions <- TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn)
  
  check1 <- nrow(predictions) == (length(xIn) * 2*(nrow(fitsIn) + sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit,na.rm = TRUE)))
  check2 <- setdiff(colnames(fitsIn), colnames(predictions)) == "fittedModel"
  check3 <- all(setdiff(colnames(predictions), colnames(fitsIn)) == c("x", "replicate", "y", "condition"))
  check4 <- all(unique(predictions$replicate) == c("Replicate1", "Replicate2"))
  check5 <- all(unique(predictions$condition) == c("null model", "Vehicle", "Treatment"))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
})

test_that(desc="allOk_H0", code={
  
  fitsIn <- splineFits %>% filter(testHypothesis == "null")
  xIn <- xNew

  predictions <- TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn)
  
  check1 <- nrow(predictions) == (length(xNew) * nrow(fitsIn))
  check2 <- setdiff(colnames(fitsIn), colnames(predictions)) == "fittedModel"
  check3 <- all(setdiff(colnames(predictions), colnames(fitsIn)) == c("x", "y"))
  check4 <- !("condition" %in% colnames(predictions))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="allOk_H1", code={
  
  fitsIn <- splineFits %>% filter(testHypothesis == "alternative")
  xIn <- xNew
  
  predictions <- TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn)
  
  check1 <- nrow(predictions) == (length(xIn) * (nrow(fitsIn) + sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit,na.rm = TRUE)))
  check2 <- setdiff(colnames(fitsIn), colnames(predictions)) == "fittedModel"
  check3 <- all(setdiff(colnames(predictions), colnames(fitsIn)) == c("x", "y", "condition"))

  expect_true(check1 & check2 & check3)
  
})

test_that(desc="modelColMissing", code={
  
  fitsIn <- splineFits %>% select(-fittedModel)
  xIn <- xNew
  
  expect_error(TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn))
  
})

test_that(desc="modelColInvalid", code={
  # If a column with assumed models is given, they are passed on to the
  # prediction. Invalid model types are handeled directly by the prediction
  # function by returning NA for each value of x.
  
  fitsIn <- splineFits %>% mutate(fittedModel = NA) # Create invalid models
  xIn <- xNew
  
  predictions <- TPP:::invoke_spline_prediction(fits = fitsIn, x = xIn)
  
  check1 <- nrow(predictions) == (length(xNew) * (nrow(fitsIn)))
  check2 <- setdiff(colnames(fitsIn), colnames(predictions)) == "fittedModel"
  check3 <- all(setdiff(colnames(predictions), colnames(fitsIn)) == c("x", "y"))
  check4 <- !("condition" %in% colnames(predictions))
  check5 <- all(is.na(predictions$y))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
})

test_that(desc="fitsMissing", code={
  
  xIn <- xNew
  
  expect_error(TPP:::invoke_spline_prediction(x = xIn))
  
})

test_that(desc="xMissing", code={
  
  fitsIn <- splineFits %>% select(-fittedModel)

  expect_error(TPP:::invoke_spline_prediction(fits = fitsIn))
  
})