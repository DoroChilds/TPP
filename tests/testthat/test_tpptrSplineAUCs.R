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

splineFits3 <- suppressMessages(
  tpptrFitSplines(data = testData, 
                  factorsH1 = c("condition", "replicate"), 
                  returnModels = TRUE, 
                  splineDF = 4,
                  nCores = 1)
)

test_that(desc="allOk", code={
  
  datIn <- testData
  fitsIn <- splineFits

  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == sum(fitsIn$testHypothesis == "null" & fitsIn$successfulFit) + 2*sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit)
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("auc", "xMin", "xMax", "condition"))
  check4 <- !all(is.na(aucs$auc))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="allOk2", code={
  
  datIn <- testData
  fitsIn <- splineFits2
  
  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == 2 * sum(fitsIn$testHypothesis == "null" & fitsIn$successfulFit) + 4*sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit)
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("replicate", "auc", "xMin", "xMax", "condition"))
  check4 <- all(unique(aucs$replicate) == c("Replicate1", "Replicate2"))
  check5 <- all(unique(aucs$condition) == c("null model", "Treatment", "Vehicle"))
  check6 <- !any(is.na(aucs$auc))
  
  expect_true(check1 & check2 & check3 & check4 & check5 & check6)
  
})

test_that(desc="allOk3", code={
  
  datIn <- testData
  fitsIn <- splineFits3
  
  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == sum(fitsIn$testHypothesis == "null" & fitsIn$successfulFit) + 4*sum(fitsIn$testHypothesis == "alternative" & fitsIn$successfulFit)
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("auc", "xMin", "xMax", "condition", "replicate"))
  check4 <- all(unique(aucs$replicate) == c("null model", "Replicate1", "Replicate2"))
  check5 <- all(unique(aucs$condition) == c("null model", "Treatment", "Vehicle"))
  # check6 <- sum(is.na(aucs$auc)) == 1 # to do: fix this so that NA dissapears
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
})

test_that(desc="allOk_H0", code={
  
  datIn <- testData
  fitsIn <- splineFits %>% filter(testHypothesis == "null")
  
  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == sum(fitsIn$successfulFit)
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("auc", "xMin", "xMax"))
  check4 <- !("condition" %in% colnames(aucs))
  # check5 <- sum(is.na(aucs$auc)) == 1 # to do: fix this so that NA dissapears
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="allOk_H1", code={
  
  datIn <- testData
  fitsIn <- splineFits %>% filter(testHypothesis == "alternative")
  
  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == sum(fitsIn$successfulFit) * 2
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("condition", "auc", "xMin", "xMax"))
  check4 <- !any(is.na(aucs$auc))
  
  expect_true(check1 & check2 & check3 & check4)

})

test_that(desc="modelColMissing", code={
  
  fitsIn <- splineFits %>% select(-fittedModel)
  datIn <- testData
  
  expect_error(tpptrSplineAUCs(data = datIn, fits = fitsIn))
  
})

test_that(desc="modelColInvalid", code={
  # If a column with assumed models is given, they are passed on to the
  # prediction. Invalid model types are handeled directly by the prediction
  # function by returning NA for each value of x.
  
  fitsIn <- splineFits %>% mutate(fittedModel = NA) # Create invalid models
  datIn <- testData
  
  aucs <- tpptrSplineAUCs(data = datIn, fits = fitsIn)
  
  check1 <- nrow(aucs) == sum(fitsIn$successfulFit)
  check2 <- setdiff(colnames(fitsIn), colnames(aucs)) == "fittedModel"
  check3 <- all(setdiff(colnames(aucs), colnames(fitsIn)) == c("auc", "xMin", "xMax"))
  check4 <- !("condition" %in% colnames(aucs))
  check5 <- all(is.na(aucs$auc))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
})