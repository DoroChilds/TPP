dataPath <- system.file("test_data", package="TPP")

dat1 <- readRDS(file.path(dataPath, "panobinostat_2D_importResults.rds"))
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_fcResults.rds"))
dat3 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults.rds")) # example input from an older experiment
dat4 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults.rds")) # example output from an older experiment (12 rows)

settings1 <- attr(dat1, "importSettings")
settings2 <- attr(dat2, "importSettings")
settings3 <- attr(dat3, "importSettings")
settings4 <- attr(dat4, "importSettings")

allCols1 <- colnames(dat1)
allCols2 <- colnames(dat2)
allCols3 <- colnames(dat3)
allCols4 <- colnames(dat4)

test_that("all_ok1", code={
  
  settingsIn <- settings1
  fieldsIn <- names(settings1) %>% setdiff(c("fcStr", "addCol", "qualColName"))
  colsIn <- colnames(dat1)
  
  colNames <- sapply(fieldsIn, function(field) {
    TPP:::checkAndReturnDataSetting(settingsIn, field, colsIn)
  }, simplify = TRUE)
  
  expect_equal(colNames, c(proteinIdCol = "representative",
                           uniqueIdCol = "unique_ID",
                           intensityStr = "sumionarea_protein_",
                           nonZeroCols = "qusm"))
  
})

test_that("all_ok2", code={
  
  settingsIn <- settings2
  fieldsIn <- names(settings2) %>% setdiff(c("addCol", "qualColName"))
  colsIn <- colnames(dat2)
  
  colNames <- sapply(fieldsIn, function(field) {
    TPP:::checkAndReturnDataSetting(settingsIn, field, colsIn)
  }, simplify = TRUE)
  
  expect_equal(colNames, c(proteinIdCol = "representative",
                           uniqueIdCol = "unique_ID",
                           intensityStr = "sumionarea_protein_",
                           nonZeroCols = "qusm",
                           fcStr = "rel_fc_protein_"))
  
})

test_that("all_ok3", code={
  
  settingsIn <- settings3
  fieldsIn <- names(settings3) %>% setdiff(c("addCol", "qualColName"))
  colsIn <- colnames(dat3)
  
  colNames <- sapply(fieldsIn, function(field) {
    TPP:::checkAndReturnDataSetting(settingsIn, field, colsIn)
  }, simplify = TRUE)
  
  expect_equal(colNames, c(proteinIdCol = "representative",
                           uniqueIdCol = "unique_ID",
                           intensityStr = "sumionarea_protein_",
                           nonZeroCols = "qusm",
                           fcStr = "rel_fc_protein_",
                           fcStrNorm = "norm_rel_fc_protein_"))
  
})

test_that("all_ok4", code={
  
  settingsIn <- settings4
  fieldsIn <- names(settings4) %>% setdiff(c("addCol", "qualColName", "slopeBounds"))
  colsIn <- colnames(dat4)
  
  colNames <- sapply(fieldsIn, function(field) {
    TPP:::checkAndReturnDataSetting(settingsIn, field, colsIn)
  }, simplify = TRUE)
  
  expect_equal(colNames, c(proteinIdCol = "representative",
                           uniqueIdCol = "Protein_ID",
                           intensityStr = "sumionarea_protein_",
                           nonZeroCols = "qusm",
                           fcStr = "rel_fc_protein_",
                           fcStrNorm = "norm_rel_fc_protein_",
                           r2Cutoff = "0.8",
                           fcCutoff = "1.5",
                           fcTolerance = "0.1"))
  
})

test_that("no_character", {
  
  settingsIn <- settings1
  settingsIn$proteinIdCol = TRUE
  fieldIn <- "proteinIdCol"
  colsIn <- colnames(dat1)
  
  expect_error(
    TPP:::checkAndReturnDataSetting(settingsIn, fieldIn, colsIn)
  )
  
})

test_that("no_numeric", {
  
  settingsIn <- settings4
  settingsIn$fcTolerance = TRUE
  fieldIn <- "fcTolerance"
  colsIn <- colnames(dat4)
  
  expect_error(
    TPP:::checkAndReturnDataSetting(settingsIn, fieldIn, colsIn)
  )
  
})


test_that("no_field", {
  
  settingsIn <- settings1
  fieldIn <- "nonsense"
  colsIn <- colnames(dat1)
  
  expect_error(
    TPP:::checkAndReturnDataSetting(settingsIn, fieldIn, colsIn)
  )
  
})

test_that("no_colname", {
  
  settingsIn <- settings1
  settingsIn$intensityStr <- "nonsense"
  fieldIn <- "intensityStr"
  colsIn <- colnames(dat1)
  
  expect_error(
    TPP:::checkAndReturnDataSetting(settingsIn, fieldIn, colsIn)
  )
  
})

