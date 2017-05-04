# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))


splineFits1 <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = c("condition"), 
                  returnModels = TRUE, splineDF = 2, nCores = 1)
)

splineFits2 <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = c("condition", "replicate"), 
                  returnModels = TRUE, splineDF = 2, nCores = 1)
)

splineFits3 <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = c("condition", "replicate"), 
                  factorsH0 = c("replicate"),
                  returnModels = TRUE, splineDF = 2, nCores = 1)
)

xNew <- seq(40, 70, by = 2)

splinePredictions1 <- TPP:::invoke_spline_prediction(fits = splineFits1, 
                                                     x = xNew)

splinePredictions2 <- TPP:::invoke_spline_prediction(fits = splineFits2, 
                                                     x = xNew)

splinePredictions3 <- TPP:::invoke_spline_prediction(fits = splineFits3, 
                                                     x = xNew)

colorBy1 <- data.frame(testHypothesis = c("alternative"), 
                       factors = c("condition"))

colorBy2 <- data.frame(testHypothesis = c("alternative", "alternative"),
                       factors = c("condition", "replicate"))

colorBy3 <- data.frame(testHypothesis = c("null", "alternative", "alternative"),
                       factors = c("condition", "condition", "replicate"))


test_that(desc="allOk1", code={
  
  datIn <- testData 
  predIn <- splinePredictions1
  fctrsIn <- colorBy1
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn,
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  check3 <- all(na.omit(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) ==
                  c( "black" , "#808080", "#da7f2d"))
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="allOk2", code={
  
  datIn <- testData 
  predIn <- splinePredictions2
  fctrsIn <- colorBy2
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn,
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  check3 <- all(na.omit(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) ==
                  c("black" , "#1B9E77", "#9B58A5", "#BBA90B", "#666666"))
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="allOk3", code={
  
  datIn <- testData 
  predIn <- splinePredictions3
  fctrsIn <- colorBy3
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn,
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  check3 <- all(na.omit(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) ==
                  c("black", "#1B9E77", "#B16548", "#D03792", "#7FA718", "#BF8B12", "#666666"))
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="sevenConditions", code={
  
  datIn <- testData %>% 
    filter(uniqueID != "CBR3") %>%
    mutate(condition = paste0(condition, as.numeric(uniqueID)),
           uniqueID = "HDAC1")
  
  fitsNew <- suppressMessages(
    tpptrFitSplines(data = datIn, factorsH1 = "condition", returnModels = TRUE, 
                    splineDF = 2, nCores = 1)
  )  
  
  fctrsIn <- colorBy1
  
  predIn <- TPP:::invoke_spline_prediction(fits = fitsNew, x = xNew)
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn, 
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- length(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) == 6
  
  expect_true(check1 & check2)
  
})

test_that(desc="onlyAlternative", code={
  
  datIn <- testData 
  
  fitsNew <- splineFits1 %>% filter(testHypothesis == "alternative")
  fctrsIn <- colorBy1
  
  predIn <- TPP:::invoke_spline_prediction(fits = fitsNew, x = xNew)
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn, 
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  check3 <- all(na.omit(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) == 
                  c("#808080", "#da7f2d"))
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="onlyNull", code={
  
  datIn <- testData 
  
  fitsNew <- splineFits3 %>% filter(testHypothesis == "null")
  fctrsIn <- colorBy3
  
  predIn <- TPP:::invoke_spline_prediction(fits = fitsNew, x = xNew)
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn, 
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  #check3 <- length(unique(ggplot2::ggplot_build(outPlot)$data[[2]]$colour)) == 2 # to do: remove the generic 'black' color
  check4 <- length(outPlot$layers) == 2
  
  expect_true(check1 & check2 & check4)
})

test_that(desc="noHypothesisNorCondition", code={
  
  datIn <- testData 
  
  fitsNew <- splineFits2 %>% filter(testHypothesis == "null") %>%
    mutate(testHypothesis = "dummy")
  
  fctrsIn <- colorBy2
  
  predIn <- TPP:::invoke_spline_prediction(fits = fitsNew, x = xNew)
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn, 
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  check3 <- length(outPlot$layers) == 2
  
  expect_true(check1 & check2 & check3)
})


test_that(desc="predictionsMissing", code={
  
  datIn <- testData
  fctrsIn <- colorBy1
  
  expect_error(
    TPP:::create_spline_plots(measurements = datIn,
                              colorBy = fctrsIn,
                              highlightIDs = c(),
                              highlightTxt = "")
  )
  
})

test_that(desc="measurementsMissing", code={
  
  predIn <- splinePredictions1
  fctrsIn <- colorBy1
  
  expect_error(
    TPP:::create_spline_plots(predictions = predIn,
                              colorBy = fctrsIn,
                              highlightIDs = c(),
                              highlightTxt = "")
  )
  
})

test_that(desc="factorsMissing", code={
  
  datIn <- testData
  predIn <- splinePredictions1
  
  expect_error(
    TPP:::create_spline_plots(measurements = datIn,
                              predictions = predIn,
                              highlightIDs = c(),
                              highlightTxt = "")
  )
  
})

test_that(desc="idColsMissingMeasurements", code={
  
  datIn <- testData %>% select(-uniqueID)
  predIn <- splinePredictions1
  fctrsIn <- colorBy1
  
  expect_error(
    TPP:::create_spline_plots(measurements = datIn, 
                              predictions = predIn,
                              colorBy = fctrsIn,
                              highlightIDs = c(),
                              highlightTxt = "")
  )
  
})

test_that(desc="idColsMissingPredictions", code={
  
  datIn <- testData 
  predIn <- splinePredictions1 %>% select(-uniqueID)
  fctrsIn <- colorBy1
  
  expect_error(
    TPP:::create_spline_plots(measurements = datIn, 
                              predictions = predIn,
                              colorBy = fctrsIn,
                              highlightIDs = c(),
                              highlightTxt = "")
  )
  
})

test_that(desc="noData", code={
  
  datIn <- testData %>% filter(uniqueID == "CBR3")
  predIn <- splinePredictions1 %>% filter(uniqueID == "CBR3")
  fctrsIn <- colorBy1
  
  outPlot <- TPP:::create_spline_plots(measurements = datIn, 
                                       predictions = predIn,
                                       colorBy = fctrsIn,
                                       highlightIDs = c(),
                                       highlightTxt = "")
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("x", "y", "colorColumn"))
  
  # Check if only the original data was added, not the failed predictions. 
  # Adding these would cause problems due to the missing color columns,
  # which is only added automatically when at least 1 successfull model fit was present:
  check3 <- length(outPlot$layers) == 1 
  
  expect_true(check1 & check2 & check3)
})
