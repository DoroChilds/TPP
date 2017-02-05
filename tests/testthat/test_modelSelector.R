testDatDummy <- expand.grid(uniqueID = paste("protein", 1:5),
                            testHypothesis = c("null", "alternative"),
                            splineDF = 3:5) %>%
  arrange(uniqueID) %>%
  mutate(aicc = 1:nrow(.))
  

test_that(desc = "allOk_defaults_minDF", code = {
  
  statsIn <- testDatDummy
  
  out <- TPP:::modelSelector(fitStats = statsIn, 
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- all(out$splineDF == 3)
  
  expect_true(check1 & check2)
  
})

test_that(desc = "allOk_defaults_maxDF", code = {
  
  statsIn <- testDatDummy %>% mutate(aicc = rev(aicc))
  
  out <- TPP:::modelSelector(fitStats = statsIn, 
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- all(out$splineDF == 5)
  
  expect_true(check1 & check2)
  
})

test_that(desc = "colsAreFactors", code = {
  
  statsIn <- testDatDummy %>% 
    mutate(uniqueID = factor(uniqueID),
           testHypothesis = factor(testHypothesis))
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- all(out$splineDF == 3)
  
  expect_true(check1 & check2)
  
})

test_that(desc = "fitStats_missing", code = {
  
  expect_error(TPP:::modelSelector(criterion = "aicc", 
                                   hypothesis = "alternative"))
  
})

test_that(desc = "criterion_not_present", code = {
  
  statsIn <- testDatDummy
  
  expect_error(TPP:::modelSelector(statsIn, criterion = "dummy",
                                   hypothesis = "alternative"))
  
})


test_that(desc = "hypothesis_col__not_present", code = {
  
  statsIn <- testDatDummy %>%select(-testHypothesis)
  
  expect_error(TPP:::modelSelector(statsIn,
                                   criterion = "aicc", 
                                   hypothesis = "alternative"))
  
})


test_that(desc = "idCol_splineDFCol_not_present", code = {
  
  statsIn <- testDatDummy %>% select(-splineDF, -uniqueID)
  
  expect_error(TPP:::modelSelector(statsIn,
                                   criterion = "aicc", 
                                   hypothesis = "alternative"))
  
})


test_that(desc = "splineDF_not_numeric", code = {
  
  statsIn <- testDatDummy %>% mutate(splineDF = as.character(splineDF))
  
  expect_error(TPP:::modelSelector(fitStats = statsIn,
                                   criterion = "aicc", 
                                   hypothesis = "alternative"))  
})


test_that(desc = "splineDF_col_allNA", code = {
  
  statsIn <- testDatDummy %>% mutate(splineDF = NA_real_)
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- all(is.na(out$splineDF))
  
  expect_true(check1 & check2)
})


test_that(desc = "splineDF_col_hasNA", code = {
  statsIn <- testDatDummy
  statsIn$splineDF[2] <- NA
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- is.na(out$splineDF[1])
  
  expect_true(check1 & check2)
})


test_that(desc = "aicc_not_numeric", code = {
  statsIn <- testDatDummy %>% mutate(aicc = as.character(aicc))
  
  expect_error(TPP:::modelSelector(fitStats = statsIn,
                                   criterion = "aicc", 
                                   hypothesis = "alternative"))    
})


test_that(desc = "aicc_ties", code = {
  
  statsIn <- testDatDummy
  statsIn$aicc[1:2] <- 1
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- out$splineDF[1] == 3
  
  expect_true(check1 & check2)
})


test_that(desc = "aicc_col_allNA", code = {
  
  statsIn <- testDatDummy
  statsIn$aicc <- NA_real_
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")   # Expected outcome: return splineDF = NA
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- all(is.na(out$splineDF))
  
  expect_true(check1 & check2)
})


test_that(desc = "aicc_col_hasNA1", code = {
  statsIn <- testDatDummy
  
  statsIn$aicc[1:6] <- NA   # Set aicc = NA for protein1
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")   # Expected outcome: return splineDF = NA only for protein1
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- is.na(out$splineDF[1])
  check3 <- all(out$splineDF[2:5] == 3)
  
  expect_true(check1 & check2 & check3)
  
})


test_that(desc = "aicc_col_hasNA2", code = {
  statsIn <- testDatDummy
  statsIn$aicc[2] <- NA
  
  out <- TPP:::modelSelector(fitStats = statsIn,
                             criterion = "aicc", 
                             hypothesis = "alternative")
  
  check1 <- nrow(out) == length(unique(statsIn$uniqueID))
  check2 <- out$splineDF[1] == 4
  check3 <- all(out$splineDF[2:5] == 3)
  
  expect_true(check1 & check2 & check3)
  
})


test_that(desc = "hypothesis_col_hasNA", code = {
  
  statsIn <- testDatDummy
  statsIn$testHypothesis[2] <- NA
  
  expect_warning(TPP:::modelSelector(fitStats = statsIn,
                                     criterion = "aicc", 
                                     hypothesis = "alternative"))
})


test_that(desc = "hypothesis_col_allNA", code = {
  statsIn <- testDatDummy
  statsIn$testHypothesis <- NA
  
  expect_error(TPP:::modelSelector(fitStats = statsIn,
                                   criterion = "aicc", 
                                   hypothesis = "alternative"))
})
