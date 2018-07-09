# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

hdacData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = hdacData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 2, nCores = 1)
)

test_that(desc="allOk", code={
  
  datIn <- hdacData
  fitsIn <- splineFits

  outPlot <- TPP:::predict_and_plot_spline_models(dat = datIn, fits = fitsIn)
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("~x", "~y", "~colorColumn"))
  
  expect_true(check1 & check2)
  
})

test_that(desc="allOk_H0", code={
  
  datIn <- hdacData
  fitsIn <- splineFits %>% filter(testHypothesis == "null")
  
  outPlot <- TPP:::predict_and_plot_spline_models(dat = datIn, fits = fitsIn)
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("~x", "~y", "~colorColumn"))
  check3 <- length(outPlot$layers) == 2 # to do: the null models should be displayed and a corresponding geom_line layer present
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="allOk_H1", code={
  
  datIn <- hdacData
  fitsIn <- splineFits %>% filter(testHypothesis != "null")
  
  outPlot <- TPP:::predict_and_plot_spline_models(dat = datIn, fits = fitsIn)
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("~x", "~y", "~colorColumn"))
  check3 <- length(outPlot$layers) == 2
  
  expect_true(check1 & check2 & check3)
  
})

test_that(desc="dataMissing", code={
  datIn <- hdacData
  fitsIn <- splineFits
  
  expect_error(TPP:::predict_and_plot_spline_models(fits = fitsIn))
  
})

test_that(desc="modelsMissing", code={
  datIn <- hdacData
  fitsIn <- splineFits
  
  expect_error(TPP:::predict_and_plot_spline_models(dat = datIn))
  
})

test_that(desc="modelColMissing", code={
  datIn <- hdacData
  fitsIn <- splineFits %>% select(-fittedModel)
  
  expect_error(TPP:::predict_and_plot_spline_models(dat = datIn, fits = fitsIn))
})

test_that(desc="modelColInvalid", code={
  # If a column with assumed models is given, they are passed on to the
  # prediction. Invalid model types are handeled directly by the prediction
  # function by returning NA for each value of x.
  
  fitsIn <- splineFits %>% mutate(fittedModel = NA) # Create invalid models
  datIn <- hdacData
  
  outPlot <- TPP:::predict_and_plot_spline_models(dat = datIn, fits = fitsIn)
  
  check1 <- inherits(outPlot, "ggplot")
  check2 <- all(paste(outPlot$mapping) == c("~x", "~y", "~colorColumn"))
  check3 <- length(outPlot$layers) == 1 # to do: the null models should be disyplayed and a corresponding geom_line layer present
  
  expect_true(check1 & check2 & check3)
  
})
