# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

hdacData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC9"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = hdacData, factorsH1 = c("condition", "replicate"), returnModels = TRUE, 
                  splineDF = 4, nCores = 1)
)

modelH0 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "null") %>%
              extract2("fittedModel"))[[1]]

modelH1 <- (splineFits %>% 
              filter(uniqueID == "HDAC1", testHypothesis == "alternative") %>%
              extract2("fittedModel"))[[1]]

test_that(desc="allOk_names_H0", code={
  
  mIn <- modelH0
  mode <- "names"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, character())
  
})

test_that(desc="allOk_values_H0", code={
  
  mIn <- modelH0
  mode <- "values"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, data.frame())
})

test_that(desc="allOk_names_H1", code={
  
  mIn <- modelH1
  mode <- "names"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, c("condition", "replicate"))
  
})

test_that(desc="allOk_values_H1", code={
  
  mIn <- modelH1
  mode <- "values"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  ref <- data.frame(condition = c("Vehicle", "Vehicle", "Treatment", "Treatment"), 
                    replicate = c("Replicate1", "Replicate2", "Replicate1", "Replicate2"),
                    stringsAsFactors = FALSE)
  
  expect_equal(fctrs, ref)
})

test_that(desc="modelMissing", code={
  
  expect_error(TPP:::extract_fit_factors(mode = "names"))
  
})

test_that(desc="modeMissing", code={
  
  mIn <- modelH1
  
  expect_error(TPP:::extract_fit_factors(splineModel = mIn))
  
})

test_that(desc="modelNULL_values", code={
  
  mIn <- NULL
  mode <- "values"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, data.frame())
  
})

test_that(desc="modelNULL_names", code={
  
  mIn <- NULL
  mode <- "names"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, character())
  
})

test_that(desc="modeWrong", code={
  
  mIn <- modelH1
  mode <- "nonsense"
  
  expect_error(TPP:::extract_fit_factors(splineModel = mIn, mode = mode))
  
})

test_that(desc="modelFitError_values", code={
  
  mIn <- try(lm(y ~ x, data = data.frame(x = NA, y = NA)), silent = TRUE)
  mode <- "values"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, data.frame())
  
})

test_that(desc="modelFitError_names", code={
  
  mIn <- try(lm(y ~ x, data = data.frame(x = NA, y = NA)), silent = TRUE)
  mode <- "names"
  
  fctrs <- TPP:::extract_fit_factors(splineModel = mIn, mode = mode)
  
  expect_equal(fctrs, character())
  
})