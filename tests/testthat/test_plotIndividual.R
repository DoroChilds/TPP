# Prepare function input:
data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))


splineFits <- suppressMessages(
  tpptrFitSplines(data = testData, factorsH1 = "condition", returnModels = TRUE, 
                  splineDF = 3:7, nCores = 1)
)

splineFits2 <- suppressMessages(
  tpptrFitSplines(data = testData, 
                  factorsH1 = c("condition", "replicate"), 
                  factorsH0 = c("replicate"),
                  returnModels = TRUE, 
                  splineDF = 3:7, nCores = 1)
)

test_that(desc="allOk_returnPlots", code={
  
  datIn <- testData
  fitsIn <- splineFits
  
  plots <- TPP:::plotIndividual(data = datIn, 
                                fittedModels = fitsIn,
                                plotAnnotation = NULL, 
                                plotDir = NULL, 
                                filePrefix  = NULL, 
                                returnPlots = TRUE, 
                                nCores = 1)
  
  check1 <- nrow(plots) == length(unique(fitsIn$uniqueID))
  check2 <- all(plots$uniqueID == unique(fitsIn$uniqueID))
  check3 <- all(sapply(plots$plot, function(p) inherits(p,"ggplot")))
  check4 <- all(colnames(plots) == c("uniqueID", "plot"))
  
  expect_true(check1 & check2 & check3 & check4)
  
})

test_that(desc="allOk_printPlots", code={
  
  datIn <- testData
  fitsIn <- splineFits
  
  paths <- TPP:::plotIndividual(data = datIn, 
                                fittedModels = fitsIn,
                                plotAnnotation = NULL, 
                                plotDir = getwd(), 
                                filePrefix = "splineFit", 
                                returnPlots = FALSE, 
                                nCores = 1)
  
  check1 <- nrow(paths) == length(unique(fitsIn$uniqueID))
  check2 <- all(paths$uniqueID == unique(fitsIn$uniqueID))
  check3 <- all(file.exists(paths$path))
  check4 <- all(colnames(paths) == c("uniqueID", "path"))
  
  unlink(file.path(getwd(), "Individual"), recursive = TRUE)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="allOk_both", code={
  
  datIn <- testData
  fitsIn <- splineFits
  
  both <- TPP:::plotIndividual(data = datIn, 
                               fittedModels = fitsIn,
                               plotAnnotation = NULL, 
                               plotDir = getwd(), 
                               filePrefix = "splineFit", 
                               returnPlots = TRUE, 
                               nCores = 1)
  
  check1 <- nrow(both) == length(unique(fitsIn$uniqueID))
  check2 <- all(both$uniqueID == unique(fitsIn$uniqueID))
  check3 <- all(sapply(both$plot, function(p) inherits(p,"ggplot")))
  check4 <- all(file.exists(both$path))
  check5 <- all(colnames(both) == c("uniqueID", "plot", "path"))
  
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
  unlink(file.path(getwd(), "Individual"), recursive = TRUE)
  
})

test_that(desc="allOk_both2", code={
  
  datIn <- testData
  fitsIn <- splineFits2
  
  both <- TPP:::plotIndividual(data = datIn, 
                               fittedModels = fitsIn,
                               plotAnnotation = NULL, 
                               plotDir = getwd(), 
                               filePrefix = "splineFit", 
                               returnPlots = TRUE, 
                               nCores = 1)
  
  check1 <- nrow(both) == length(unique(fitsIn$uniqueID))
  check2 <- all(both$uniqueID == unique(fitsIn$uniqueID))
  check3 <- all(sapply(both$plot, function(p) inherits(p,"ggplot")))
  check4 <- all(file.exists(both$path))
  check5 <- all(colnames(both) == c("uniqueID", "plot", "path"))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
  unlink(file.path(getwd(), "Individual"), recursive = TRUE)
  
})

test_that(desc="noData", code={
  
  datIn <- testData %>% filter(uniqueID == "CBR3")
  fitsIn <- splineFits
  
  plots <- TPP:::plotIndividual(data = datIn, 
                                fittedModels = fitsIn,
                                plotAnnotation = NULL, 
                                plotDir = NULL, 
                                filePrefix = NULL, 
                                returnPlots = TRUE, 
                                nCores = 1)
  p <- plots$plot[[1]]
  
  check1 <- nrow(plots) == 1
  check2 <- inherits(p, "ggplot")
  
  # Check if only the original data was added, not the failed predictions. 
  # Adding these would cause problems due to the missing color column 'condition'
  # which is only added automatically when at least 1 successfull model fit was present:
  check3 <- length(p$layers) == 1 
  
  expect_true(check1 & check2 & check3)
})