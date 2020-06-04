data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrData %>% lapply(function(d) {
  d[featureNames(d) %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"),]
})

testDataTidy <- tpptrTidyUpESets(testData, returnType = "exprs")

resPath <- getwd()

expectedCols <- c("F_statistic", "F_moderated",  "F_scaled", "residual_df_H1", 
                  "prior_df_H1", "df1", "df2", "df2_moderated", "posterior_var_H1",
                  "p_NPARC", "p_adj_NPARC")

test_that(desc="allOk1", code={
  
  datIn <- testDataTidy
  facH0 <- character()
  facH1 <- c("condition")
  doPlot <- FALSE
  cores <- 1
  dfs <- 3
  addCol <- NULL
  
  out <- tpptrSplineFitAndTest(data = datIn,
                               factorsH1 = facH1,
                               factorsH0 = facH0,
                               resultPath = resPath,
                               doPlot = doPlot,
                               nCores = cores, 
                               splineDF = dfs, 
                               additionalCols = addCol)
  
  check1 <- all(unique(out$Protein_ID) == unique(datIn$uniqueID))
  check2 <- all(expectedCols %in% colnames(out))
  check3 <- !all(apply(select(out,!!expectedCols), 2, is.na))
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="allOk2", code={
  
  datIn <- testDataTidy
  facH0 <- character()
  facH1 <- c("condition", "replicate")
  doPlot <- FALSE
  cores <- 1
  dfs <- 3
  addCol <- NULL
  
  out <- tpptrSplineFitAndTest(data = datIn,
                               factorsH1 = facH1,
                               factorsH0 = facH0,
                               resultPath = resPath,
                               doPlot = doPlot,
                               nCores = cores, 
                               splineDF = dfs, 
                               additionalCols = addCol)
  
  check1 <- all(unique(out$Protein_ID) == unique(datIn$uniqueID))
  check2 <- all(expectedCols %in% colnames(out))
  check3 <- !all(apply(select(out,!!expectedCols), 2, is.na))
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="allOk3", code={
  
  datIn <- testDataTidy
  facH0 <- c("replicate")
  facH1 <- c("condition", "replicate")
  doPlot <- FALSE
  cores <- 1
  dfs <- 3
  addCol <- NULL
  
  out <- tpptrSplineFitAndTest(data = datIn,
                               factorsH1 = facH1,
                               factorsH0 = facH0,
                               resultPath = resPath,
                               doPlot = doPlot,
                               nCores = cores, 
                               splineDF = dfs, 
                               additionalCols = addCol)
  
  check1 <- all(unique(out$Protein_ID) == unique(datIn$uniqueID))
  check2 <- all(expectedCols %in% colnames(out))
  check3 <- !all(apply(select(out,!!expectedCols), 2, is.na))
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="allOk3_doplot", code={
  
  datIn <- testDataTidy
  facH0 <- c("replicate")
  facH1 <- c("condition", "replicate")
  doPlot <- TRUE
  cores <- 1
  dfs <- 3
  addCol <- NULL
  
  out <- tpptrSplineFitAndTest(data = datIn,
                               factorsH1 = facH1,
                               factorsH0 = facH0,
                               resultPath = resPath,
                               doPlot = doPlot,
                               nCores = cores, 
                               splineDF = dfs, 
                               additionalCols = addCol)
  
  check1 <- all(unique(out$Protein_ID) == unique(datIn$uniqueID))
  check2 <- all(expectedCols %in% colnames(out))
  check3 <- !all(apply(select(out,!!expectedCols), 2, is.na))
  check4 <- out$splinefit_plot %>% file.path(resPath, .) %>% file.exists %>% all
  
  unlink(file.path(resPath, "Spline_Fits"), recursive = TRUE)
  unlink(file.path(resPath, "QCplots_fTest.pdf"), recursive = TRUE)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="allOk_eSets", code={
  
  datIn <- testData
  facH0 <- character()
  facH1 <- c("condition")
  doPlot <- FALSE
  cores <- 1
  dfs <- 3
  addCol <- NULL
  
  out <- tpptrSplineFitAndTest(data = datIn,
                               factorsH1 = facH1,
                               factorsH0 = facH0,
                               resultPath = resPath,
                               doPlot = doPlot,
                               nCores = cores, 
                               splineDF = dfs, 
                               additionalCols = addCol)
  
  check1 <- all(unique(out$Protein_ID) == unique(datIn$uniqueID))
  check2 <- all(expectedCols %in% colnames(out))
  check3 <- !all(apply(select(out,!!expectedCols), 2, is.na))
  
  expect_true(check1 & check2 & check3)
})