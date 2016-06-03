dirExmpDat <- system.file("example_data", package="TPP")

data(hdacTR_smallExample)
data(hdacCCR_smallExample)
refLabels <- c('126', '127L', '127H', '128L', '128H','129L', '129H','130L', '130H', '131L')
refTemps  <- c(67, 63, 59, 56, 53, 50, 47, 44, 41, 37)
refTempMat <- matrix(rep(refTemps,4), nrow=4, byrow=TRUE, dimnames = list(c(1:4), refLabels))
refNames <- c("Vehicle_1", "Vehicle_2", "Panobinostat_1", "Panobinostat_2")
refConds <- c("Vehicle", "Vehicle", "Treatment", "Treatment")
refComps <- c(ComparisonVT1="Panobinostat_1_vs_Vehicle_1",
              ComparisonVT2="Panobinostat_2_vs_Vehicle_2")

## ------------------------------------------------------------------------- ##
## function 'importCheckConfigTable':
## ------------------------------------------------------------------------- ##
test_that(desc="confgCheck", code={
  message("Checking function importCheckConfigTable")
  confgList <- importCheckConfigTable(infoTable=hdacTR_config, type="TR")
  refList   <- list(expNames = refNames,
                    expCond  = refConds,
                    files    = NULL,
                    compStrs = refComps,
                    labels   = refLabels,
                    tempMatrix = refTempMat)
  expect_equal(confgList, refList)
})

test_that(desc="confgCheck_expColNULL", code={
  hdacTR_config$Experiment <- NULL
  expect_error(importCheckConfigTable(infoTable=hdacTR_config, type="TR"))
})

test_that(desc="'confgCheck_expColConvertAlnum", code={
  oldCol <- paste(hdacTR_config$Experiment, c("", ":", "'", "!"))
  hdacTR_config$Experiment <- oldCol
  confgList <- importCheckConfigTable(infoTable=hdacTR_config, type="TR")
  newCol <- confgList$expNames
  expect_equal(newCol, gsub("([^[:alnum:]])", "_", oldCol))
})

test_that(desc="'confgCheck_expColEmptyEntries", code={
  hdacTR_config$Experiment <- gsub("Vehicle_1", "", hdacTR_config$Experiment)
  confgList <- importCheckConfigTable(infoTable=hdacTR_config, type="TR")
  expect_equal(confgList$expNames, hdacTR_config$Experiment[-1])
})




## ------------------------------------------------------------------------- ##
## function 'importFct_readConfigTable':
## ------------------------------------------------------------------------- ##
test_that(desc='confgImportDF', code={
  message("Checking function importFct_readConfigTable (import TR-config from data frame).")
  cfg <- importFct_readConfigTable(cfg=hdacTR_config)
  expect_equal(cfg, hdacTR_config)
})

# test_that(desc='confgImportFileTR', code={
#   message("Checking function importFct_readConfigTable (import TR-config from xlsx).")
#   cfg <- importFct_readConfigTable(cfg=file.path(dirExmpDat, "TR_example_data", "Panobinostat_TPP-TR_config.xlsx"))
#   cfg$Path <- NULL
#   expect_equal(cfg, hdacTR_config)
# })

test_that(desc='confgImportFileCCR', code={
  message("Checking function importFct_readConfigTable (import CCR-config from xlsx).")
  cfg <- importFct_readConfigTable(cfg=file.path(dirExmpDat, "CCR_example_data", "Panobinostat_TPP-CCR_config.xlsx"))
  cfg$Path <- NULL
  expect_equal(cfg, hdacCCR_config)
})


## ------------------------------------------------------------------------- ##
## function 'importFct_checkExperimentCol':
## ------------------------------------------------------------------------- ##
test_that(desc="'confgExpCol_allok", code={
  oldCol <- c("abc", "def", "123")
  newCol <- importFct_checkExperimentCol(expCol=oldCol)
  expect_equal(newCol, oldCol)
})

test_that(desc="'confgExpCol_convertAlnum", code={
  oldCol <- c("ok1", "abc!", "def'", "123$", "ok2")
  newCol <- importFct_checkExperimentCol(expCol=oldCol)
  expect_equal(newCol, gsub("([^[:alnum:]])", "_", oldCol))
})

test_that(desc="'confgExpCol_NULLcol", code={
  oldCol <- NULL
  expect_error(importFct_checkExperimentCol(expCol=oldCol))
})


## ------------------------------------------------------------------------- ##
## function 'importFct_checkConditions':
## ------------------------------------------------------------------------- ##
## test 1: correct condition input
test_that(desc='checkCondCorrect', code={
  message("Checking function importFct_checkConditions.")
  newConds <- importFct_checkConditions(refConds)
  expect_equal(newConds, refConds)
})

## test 2: NULL input -> create NAs, produce message
test_that(desc='checkCondNull', code={
  message("Checking function importFct_checkConditions.")
  newConds <- importFct_checkConditions(NULL, 4)
  expect_equal(newConds, rep(NA_character_, 4))
})

## test 3: Conditions named differently than 'Vehicle' and 'Treatment' -> create NAs, produce message
test_that(desc='checkCond_wrongChars', code={
  message("Checking function importFct_checkConditions.")
  newConds <- importFct_checkConditions(c("A","B", "C", "D"), 4)
  expect_equal(newConds, rep(NA_character_, 4))
})

## test 4: Condition column contains numbers instead of characters -> create NAs, produce message
test_that(desc='checkCond_noChars', code={
  message("Checking function importFct_checkConditions.")
  newConds <- importFct_checkConditions(1:4, 4)
  expect_equal(newConds, rep(NA_character_, 4))
})

## ------------------------------------------------------------------------- ##
## function 'importFct_checkComparisons':
## ------------------------------------------------------------------------- ##
test_that(desc = 'checkComps_correct', code = {
  message("Checking function importFct_checkComparisons.")
  compStrs <- importFct_checkComparisons(confgTable = hdacTR_config)
  refStrs <- c(ComparisonVT1 = "Panobinostat_1_vs_Vehicle_1",
               ComparisonVT2 = "Panobinostat_2_vs_Vehicle_2")
  expect_equal(compStrs, refStrs)
})

test_that(desc = 'checkComps_nullCol', code = {
  message("Checking function importFct_checkComparisons.")
  testtab <- hdacTR_config
  testtab$ComparisonVT1 <- testtab$ComparisonVT2 <- NULL
  compStrs <- importFct_checkComparisons(confgTable = testtab)
  expect_equal(compStrs, NULL)
})

## Still to check:
## ------------------------------------------------------------------------- ##
## function 'importCheckTemperatures':
## -------------------------------------------------------------------------- ##
test_that(desc = 'checkTemperatures_all_ok', code = {
  message("Checking function checkTemperatures")
  ## Input:
  T <- hdacTR_config[,refLabels]
  ## Expected output:
  T_ref <- as.matrix(T)
  ## Test:
  T_new <- importCheckTemperatures(temp = T)
  expect_equal(T_ref, T_new)  
})

test_that(desc = 'checkTemperatures_NA_values', code = {
  message("Checking function checkTemperatures")
  ## Input:
  T <- hdacTR_config[,refLabels]
  T[1:2,1:5] <- NA
  ## Expected output:
  T_ref <- as.matrix(T)
  ## Test:
  T_new <- importCheckTemperatures(temp = T)
  expect_equal(T_ref, T_new)  
})

test_that(desc = 'checkTemperatures_NA_rows', code = {
  message("Checking function checkTemperatures")
  ## Input:
  T <- hdacTR_config[,refLabels]
  T[1:2,] <- NA
  ## Test:
  expect_error(importCheckTemperatures(temp = T))
})

# importCheckExperimentNames

# importCheckDataFormat

## ------------------------------------------------------------------------- ##
## function 'importFct_makeOutputDirs':
## -------------------------------------------------------------------------- ##
test_that(desc = 'checkOutDirs_path_and_names_given', code = {
  message("Checking function importFct_makeOutputDirs")
  ## Input:
  d <- getwd(); 
  f <- file.path(d, "dummy.txt")
  ## Expected output:
  e <- list(doWrite=TRUE, pathDataObj=file.path(d, "dataObj"), outDir=d)
  ## Test:
  res <- importFct_makeOutputDirs(outDir=d, fNames=f)
  unlink(res$pathDataObj, recursive=TRUE)
  expect_equal(res, e)
})

test_that(desc = 'checkOutDirs_path_NULL', code = {
  message("Checking function importFct_makeOutputDirs")
  ## Input:
  f <- file.path(getwd(), "dummy.txt")
  ## Expected output:
  newd <- file.path(getwd(), "TPP_results")
  e <- list(doWrite=TRUE, pathDataObj=file.path(newd, "dataObj"), outDir=newd)
  ## Test:
  res <- importFct_makeOutputDirs(outDir=NULL, fNames=f)
  unlink(res$outDir, recursive=TRUE)
  expect_equal(res, e)
})

test_that(desc = 'checkOutDirs_files_NULL', code = {
  message("Checking function importFct_makeOutputDirs")
  ## Input:
  d <- getwd()
  ## Expected output:
  e <- list(doWrite=TRUE, pathDataObj=file.path(d, "dataObj"), outDir=d)
  ## Test:
  res <- importFct_makeOutputDirs(outDir=d, fNames=NULL)
  unlink(res$pathDataObj, recursive=TRUE)
  expect_equal(res, e)
})

test_that(desc = 'checkOutDirs_both_NULL', code = {
  message("Checking function importFct_makeOutputDirs")
  ## Expected output:
  e <- list(doWrite=FALSE, pathDataObj=NULL, outDir=NULL)
  ## Test:
  res <- importFct_makeOutputDirs(outDir=NULL, fNames=NULL)
  expect_equal(res, e)
})