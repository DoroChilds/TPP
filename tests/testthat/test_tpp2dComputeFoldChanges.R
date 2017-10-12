dataPath <- system.file("test_data", package="TPP")

# Load function input:
dat1 <- readRDS(file.path(dataPath, "panobinostat_2D_importResults.rds")) # example input from an older experiment
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_importResults2.rds")) # example input from an older experiment

# Load expected result for the given input:
out1 <- readRDS(file.path(dataPath, "panobinostat_2D_fcResults.rds"))
out2 <- readRDS(file.path(dataPath, "panobinostat_2D_fcResults2.rds"))

# Start tests:
test_that(desc="all_ok1", code={
  
  datIn <- dat1
  fcStrIn <- "rel_fc_protein_"
  
  new1 <- tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStrIn)
  
  expect_identical(out1, new1)
})

test_that(desc="all_ok2", code={
  
  datIn <- dat2
  fcStrIn <- "rel_fc_protein_"
  
  new2 <- tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStrIn)
  
  expect_identical(out2, new2)
})

test_that(desc="default_fcStr", code={
  
  datIn <- dat1
  
  new1 <- tpp2dComputeFoldChanges(data = datIn)
  
  new_fc_cols <- c("rel_fc_5", "rel_fc_1", "rel_fc_0.143", "rel_fc_0.02", 
                   "rel_fc_0")
  
  ref_fc_cols <- gsub("rel_fc", "rel_fc_protein", new_fc_cols)
  
  check1 <- all(new_fc_cols %in% colnames(new1))
  check2 <- all(new1[, new_fc_cols] == out1[, ref_fc_cols], na.rm = TRUE)
  
  expect_true(check1 & check2)
})


test_that(desc="no_intensities", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$intensityStr <- NULL
  
  fcStrIn <- "rel_fc_protein_"
  
  expect_error(
    tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStrIn)
  )  
})

test_that(desc="intensities_not_character", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$intensityStr <- numeric()
  
  fcStrIn <- "rel_fc_protein_"
  
  expect_error(
    tpp2dComputeFoldChanges(data = datIn, newFcStr = fcStrIn)
  )  
})

test_that(desc="intensities_not_found", code={
  
  datIn <- dat1
  attr(datIn, "importSettings")$intensityStr <- "dummy"
  
  expect_error(
    tpp2dComputeFoldChanges(data = datIn)
  )  
})

test_that(desc="warning_deprecated_fct_arg1", code={
  
  datIn <- dat1

  expect_warning(
    tpp2dComputeFoldChanges(data = datIn,
                            intensityStr = "sumion_area")
  )
})

test_that(desc="warning_deprecated_fct_arg2", code={
  
  datIn <- dat1

  expect_warning(
    tpp2dComputeFoldChanges(data = datIn,
                            fcStr = fcStrIn)
  )
})

test_that(desc="warning_deprecated_fct_arg3", code={
  
  datIn <- dat1

  expect_warning(
    tpp2dComputeFoldChanges(configTable = NULL, data = datIn)
  )
})

test_that(desc="data_missing", code={
  
  expect_error(
    tpp2dComputeFoldChanges()
  )
})


