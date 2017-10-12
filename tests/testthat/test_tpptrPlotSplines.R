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

testStats <- tpptrFTest(fittedModels = splineFits)

testStats2 <- tpptrFTest(fittedModels = splineFits2)

resPath <- getwd()

test_that(desc="allOk", code={
  
  datIn <- testData
  fitsIn <- splineFits
  testsIn <- testStats
  
  paths <- tpptrPlotSplines(data = datIn, 
                            fittedModels = fitsIn, 
                            testResults = testsIn, 
                            resultPath = resPath,
                            individual = TRUE,
                            control = list(nCores = 1))
  
  check1 <- all(names(paths) == c("individual", "overview"))
  check2 <- nrow(paths$individual) == length(unique(fitsIn$uniqueID))
  check3 <- all(paths$individual$uniqueID == unique(fitsIn$uniqueID))
  check4 <- paths$individual$path %>% file.path(resPath, .) %>% file.exists %>% all
  check5 <- all(colnames(paths$individual) == c("uniqueID", "path"))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
  unlink(file.path(resPath, "Spline_Fits"), recursive = TRUE)
})

test_that(desc="allOk2", code={
  
  datIn <- testData
  fitsIn <- splineFits2
  testsIn <- testStats2
  
  paths <- tpptrPlotSplines(data = datIn, 
                            fittedModels = fitsIn, 
                            testResults = testsIn, 
                            resultPath = resPath,
                            individual = TRUE,
                            control = list(nCores = 1))
  
  check1 <- all(names(paths) == c("individual", "overview"))
  check2 <- nrow(paths$individual) == length(unique(fitsIn$uniqueID))
  check3 <- all(paths$individual$uniqueID == unique(fitsIn$uniqueID))
  check4 <- paths$individual$path %>% file.path(resPath, .) %>% file.exists %>% all
  check5 <- all(colnames(paths$individual) == c("uniqueID", "path"))
  
  expect_true(check1 & check2 & check3 & check4 & check5)
  
  unlink(file.path(resPath, "Spline_Fits"), recursive = TRUE)
})

test_that(desc="dataMissing", code={
  
  fitsIn <- splineFits
  testsIn <- testStats
  resPath <- getwd()
  
  expect_error(
    tpptrPlotSplines(fittedModels = fitsIn, 
                     testResults = testsIn, 
                     resultPath = resPath)
  )
})

test_that(desc="modelsMissing", code={
  
  datIn <- testData
  testsIn <- testStats
  resPath <- getwd()
  
  expect_error(
    tpptrPlotSplines(data = datIn,
                     testResults = testsIn, 
                     resultPath = resPath)
  )
})

test_that(desc="testStatsMissing", code={
  
  datIn <- testData
  fitsIn <- splineFits
  resPath <- getwd()
  
  expect_error(
    tpptrPlotSplines(data = datIn,
                     fittedModels = fitsIn, 
                     resultPath = resPath)
  )
})

test_that(desc="pathMissing", code={
  
  datIn <- testData
  fitsIn <- splineFits
  testsIn <- testStats
  
  expect_warning(
    tpptrPlotSplines(data = datIn,
                     fittedModels = fitsIn,
                     testResults = testsIn,
                     control = list(nCores = 1))
  )
})

test_that(desc="argDeprecated1", code={
  datIn <- testData
  fitsIn <- splineFits
  testsIn <- testStats
  resPath <- getwd()
  
  expect_warning({
    paths <- tpptrPlotSplines(data = datIn, 
                              factorsH1 = "condition",
                              fittedModels = fitsIn, 
                              testResults = testsIn, 
                              resultPath = resPath,
                              control = list(nCores = 1))
    
    unlink(file.path(resPath, "Spline_Fits"), recursive = TRUE)
  })
})

test_that(desc="argDeprecated2", code={
  datIn <- testData
  fitsIn <- splineFits
  testsIn <- testStats
  resPath <- getwd()
  
  expect_warning({    
    paths <- tpptrPlotSplines(data = datIn, 
                              factorsH0 = "condition",
                              fittedModels = fitsIn, 
                              testResults = testsIn, 
                              resultPath = resPath,
                              control = list(nCores = 1))
    
    unlink(file.path(resPath, "Spline_Fits"), recursive = TRUE)
  })
})