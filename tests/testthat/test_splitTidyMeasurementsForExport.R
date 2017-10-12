# Currently only tests NPARC scenario
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

test_that("allOK", {
  
  datIn <- testData
  
  out <- TPP:::splitTidyMeasurementsForExport(measurements = datIn, 
                                              proteinInfos = NULL) 
  
  expectedCols <- c("fcDF", "curveParDF", "plotCol", "presenceDF", "modelInfoDF", "otherAnnotDF")
  
  check1 <- all(names(out) == expectedCols)
  check2 <- "Protein_ID" %in% colnames(out$fcDF)
  check3 <- nrow(out$fcDF) == length(unique(datIn$uniqueID))
  
  expect_true(check1 & check2 & check3)
})

# to do: add more tests here