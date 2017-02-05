# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

hdacData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC9"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = hdacData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 4, nCores = 1)
)

xNew <- 1:100

modelH0 <- (splineFits %>% 
               filter(uniqueID == "HDAC1", testHypothesis == "null") %>%
               extract2("fittedModel"))[[1]]

modelH1 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "alternative") %>%
              extract2("fittedModel"))[[1]]

test_that(desc="allOk_H0", code={
  
  mIn <- modelH0
  xIn <- xNew
  
  prediction <- TPP:::predict_spline(splineModel = mIn, x = xIn)
  
  check1 <- nrow(prediction) == length(xNew)
  check2 <- all(colnames(prediction) == c("x", "y"))
  check3 <- all(prediction$x == xNew)
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="allOk_H1", code={
  
  mIn <- modelH1
  xIn <- xNew
  
  prediction <- TPP:::predict_spline(splineModel = mIn, x = xIn)
  
  check1 <- nrow(prediction) == (2 * length(xNew))
  check2 <- all(colnames(prediction) == c("x", "condition", "y"))
  check3 <- all(prediction$x == rep(xNew, each = 2))
  check4 <- all(unique(prediction$condition) == c("Vehicle", "Treatment"))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="modelMissing", code={
  
  expect_error(TPP:::predict_spline(x = 1:10))
  
})

test_that(desc="xMissing", code={
  
  mIn <- modelH0
  
  expect_error(TPP:::predict_spline(splineModel = mIn))
})

test_that(desc="modelNULL", code={
  
  mIn <- NULL
  xIn <- xNew
  
  prediction <- TPP:::predict_spline(splineModel = mIn, x = xIn)
  
  check1 <- nrow(prediction) == length(xNew)
  check2 <- all(colnames(prediction) == c("x", "y"))
  check3 <- all(prediction$x == xNew)
  check4 <- all(is.na(prediction$y))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="xNULL", code={
  mIn <- modelH0
  xIn <- c() # equivalent to NULL
  
  expect_error(TPP:::predict_spline(splineModel = mIn, x = xIn))
  
})

test_that(desc="modelFitError", code={
  
  mIn <- try(lm(y ~ x, data = data.frame(x = NA, y = NA)), silent = TRUE)
  xIn <- xNew
  
  prediction <- TPP:::predict_spline(splineModel = mIn, x = xIn)
  
  check1 <- nrow(prediction) == length(xNew)
  check2 <- all(colnames(prediction) == c("x", "y"))
  check3 <- all(prediction$x == xNew)
  check4 <- all(is.na(prediction$y))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="xNonNumeric", code={
  
  mIn <- modelH0
  xIn <- as.character(xNew)
  
  expect_error(TPP:::predict_spline(splineModel = mIn, x = xIn))
})

test_that(desc="xEmpty", code={
  mIn <- modelH0
  xIn <- numeric()
  
  expect_error(TPP:::predict_spline(splineModel = mIn, x = xIn))
  
})

