# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 3:4, nCores = 1)
)

splineFits2 <- suppressMessages(
  tpptrFitSplines(data = testData, 
                  factorsH1 = c("condition", "replicate"), 
                  factorsH0 = c("replicate"),
                  returnModels = TRUE, 
                  splineDF = 3:4,
                  nCores = 1)
)

modelH0 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "null") %>%
              extract2("fittedModel"))[[1]]

modelH1_1 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "alternative") %>%
              extract2("fittedModel"))[[1]]

modelH1_2 <- (splineFits2 %>% 
                filter(uniqueID == "HDAC1", testHypothesis == "alternative") %>%
                extract2("fittedModel"))[[1]]

test_that(desc="allOk_H0", code={
  
  mIn <- modelH0
  xmin <- 30
  xmax <- 60
  out <- TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax)
  
  check1 <- is.data.frame(out)
  check2 <- nrow(out) == 1
  check3 <- all(colnames(out) == c("auc"))
  check4 <- all(out$auc > 0)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="allOk_H1_1", code={
  
  mIn <- modelH1_1
  xmin <- 30
  xmax <- 60
  out <- TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax)
  
  check1 <- is.data.frame(out)
  check2 <- nrow(out) == 2
  check3 <- all(colnames(out) == c("condition", "auc"))
  check4 <- all(out$auc > 0)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="allOk_H1_2", code={
  
  mIn <- modelH1_2
  xmin <- 30
  xmax <- 60
  out <- TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax)
  
  check1 <- is.data.frame(out)
  check2 <- nrow(out) == 4
  check3 <- all(colnames(out) == c("condition", "replicate", "auc"))
  check4 <- all(out$auc > 0)
  
  expect_true(check1 & check2 & check3 & check4)
})


test_that(desc="modelMissing", code={
  
  xmin <- 30
  xmax <- 60
  expect_error(TPP:::compute_spline_auc(xmin = xmin, xmax = xmax))
  
})

test_that(desc="xMinMissing", code={
  
  mIn <- modelH0
  xmax <- 60
  
  expect_error(TPP:::compute_spline_auc(splineModel = mIn, xmax = xmax))
})

test_that(desc="xMaxMissing", code={
  
  mIn <- modelH0
  xmin <- 30

  expect_error(TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin))
})

test_that(desc="modelNULL", code={
  
  mIn <- NULL
  xmin <- 30
  xmax <- 60
  
  out <- TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax)
  
  check1 <- is.data.frame(out)
  check2 <- nrow(out) == 1
  check3 <- all(colnames(out) == c("auc"))
  check4 <- is.na(out$auc)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="xNULL", code={
  
  mIn <- modelH0
  xmin <- c() # equivalent to NULL
  xmax <- 60
  
  expect_error(TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax))
  
})

test_that(desc="modelFitError", code={
  
  mIn <- try(lm(y ~ x, data = data.frame(x = NA, y = NA)), silent = TRUE)
  xmin <- 30
  xmax <- 60
  
  out <- TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax)
  
  check1 <- is.data.frame(out)
  check2 <- nrow(out) == 1
  check3 <- all(colnames(out) == c("auc"))
  check4 <- is.na(out$auc)
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="xNonNumeric", code={
  mIn <- modelH0
  xmin <- "30"
  xmax <- 60
  expect_error(TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax))
})

test_that(desc="xEmpty", code={
  mIn <- modelH0
  xmin <- numeric()
  xmax <- 60
  
  expect_error(TPP:::compute_spline_auc(splineModel = mIn, xmin = xmin, xmax = xmax))
  
})