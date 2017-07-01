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

splineStats <- tpptrFTest(fittedModels = splineFits) %>% rename(Protein_ID = uniqueID)

dataList <- TPP:::splitTidyMeasurementsForExport(measurements = testData, 
                                                 proteinInfos = NULL)


test_that(desc = "allOk", code = {
  
  datIn <- dataList
  statsIn <- splineStats
  
  out <- TPP:::mergeOutputTables_TR(dataList = datIn, pValDF = statsIn, qualCheckDF = NULL)
  
  check1 <- nrow(out) == length(dataList$fcDF$Protein_ID)
  check2 <- ncol(out) == 52
  check3 <- identical(out$p_NPARC, statsIn$p_NPARC)
  
  expect_true(check1 & check2 & check3)
})

test_that(desc = "dataMissing", code = {
  
  statsIn <- splineStats
  
  expect_error(TPP:::mergeOutputTables_TR(pValDF = statsIn, qualCheckDF = NULL))
})

test_that(desc = "fieldMissing", code = {
  
  datIn <- dataList %>% inset2("fcDF", NULL)
  statsIn <- splineStats
  
  expect_error(TPP:::mergeOutputTables_TR(dataList = datIn, pValDF = statsIn, qualCheckDF = NULL))
})