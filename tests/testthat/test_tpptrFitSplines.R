data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

testData <- tpptrTidyUpESets(tpptrData, returnType = "exprs") %>%
  filter(uniqueID %in% c("HDAC1", "HDAC2", "HDAC9", "CBR3"))


test_that(desc="allOk1", code={
  
  datIn <- testData
  facH0 <- character(0)
  facH1 <- c("condition")
  cores <- 1
  dfs <- 3:4
  
  out <- tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                         splineDF = dfs, returnModels = FALSE, nCores = cores)
  
  check1 <- all(unique(out$uniqueID) == unique(datIn$uniqueID))
  check2 <- sum(out$successfulFit, na.rm = TRUE) == 5
  
  expect_true(check1 & check2)
  
})

test_that(desc="allOk2", code={
  
  datIn <- testData
  facH0 <- c("replicate")
  facH1 <- c("condition", "replicate")
  cores <- 1
  dfs <- 3:4
  
  out <- tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                         splineDF = dfs, returnModels = FALSE, nCores = cores)
  
  check1 <- all(unique(out$uniqueID) == unique(datIn$uniqueID))
  check2 <- sum(out$successfulFit, na.rm = TRUE) == 4
  
  expect_true(check1 & check2)
})

test_that(desc="allOK_AUC", code={
  
  datIn <- testData
  facH0 <- character(0)
  facH1 <- c("condition")
  cores <- 1
  dfs <- 3:4
  
  expect_warning(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                                 splineDF = dfs, computeAUC = TRUE, 
                                 returnModels = FALSE, nCores = cores))
  
})

test_that(desc="allOk_modelReturn", code={
  
  datIn <- testData
  facH0 <- character(0)
  facH1 <- c("condition")
  cores <- 1
  dfs <- 3:4
  
  out <- tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                         splineDF = dfs, returnModels = TRUE, nCores = cores)
  
  check1 <- all(unique(out$uniqueID) == unique(datIn$uniqueID))
  check2 <- sum(out$successfulFit, na.rm = TRUE) == 5
  check3 <- "fittedModel" %in% colnames(out)
  check4 <- sum(out$fittedModel %>% sapply(function(l) class(l)) == "lm", na.rm = TRUE) == 4
  
  expect_true(check1 & check2 & check3 & check4)
  
})

# test_that(desc="allOk_Parallel", code={
#   
#   datIn <- tpptrTidyUpESets(tpptrData, returnType = "exprs")
#   facH0 <- character(0)
#   facH1 <- c("condition")
#   dfs <- 3:6
#   
#   t1 <- system.time(
#     out <- tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
#                            splineDF = dfs, returnModels = TRUE, nCores = 4)
#   )
#   
#   t2 <- system.time(
#     out <- tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
#                            splineDF = dfs, returnModels = TRUE, nCores = 1)
#   )
#   
#   expect_true(t1["elapsed"] < t2["elapsed"])
#   
# })

test_that(desc="dataMissing", code={
  
  facH0 <- c()
  facH1 <- c("condition")
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
  
})

test_that(desc="factorH1NonCharacter", code={
  
  datIn <- testData
  facH1 <- factor(c("condition"))
  facH0 <- character()
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="factorH0NonCharacter", code={
  
  datIn <- testData
  facH1 <- "condition"
  facH0 <- factor(0)
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="dfsNonNumeric", code={
  
  datIn <- testData
  facH1 <- "condition"
  facH0 <- character(0)
  cores <- 1
  dfs <- as.character(3:4)
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="factorH1Missing", code={
  
  datIn <- testData
  facH0 <- character()
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="factorH1NotInData", code={
  
  datIn <- testData
  facH1 <- c("condition", "nonsense")
  facH0 <- character()
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="factorH1AllNA", code={
  
  datIn <- testData %>% mutate(nonsense = NA)
  facH1 <- c("condition", "replicate", "nonsense")
  facH0 <- character()
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})

test_that(desc="factorH1only1Level", code={
  
  datIn <- testData %>% mutate(condition = 1, replicate = 2)
  facH1 <- c("condition", "replicate")
  facH0 <- character()
  cores <- 1
  dfs <- 3:4
  
  expect_error(tpptrFitSplines(data = datIn, factorsH1 = facH1, factorsH0 = facH0, 
                               splineDF = dfs, returnModels = FALSE, nCores = cores))
  
})
