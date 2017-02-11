dataPath <- system.file("test_data", package="TPP")

# Load function input:
dat1 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults.rds")) # example input
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults2.rds")) # example input from an older experiment (12 rows)
dat3 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults3.rds")) # example input from an older experiment (20 rows)

# Load expected result for the given input:
out1 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults.rds")) # example output
out2 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults2.rds")) # example output from an older experiment (12 rows)
out3 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults3.rds")) # example output from an older experiment (20 rows)

# Start tests:
test_that("all_ok1", code={
  datIn <- dat1
  
  new <- tpp2dCurveFit(data = datIn, nCores = 1) %>% select(-row_derived_from_non_unique_identifiers)

  expect_equal(new %>% mutate(pEC50_quality_check = as.numeric(pEC50_quality_check)), 
               out1 %>% mutate(pEC50_quality_check = as.numeric(pEC50_quality_check)), 
               tolerance = 1e-5) # There are some small deviations in the estimated slopes on linux, which do not occur on the Mac.
})

test_that("all_ok2", code={
  datIn <- dat2
  
  new <- tpp2dCurveFit(data = datIn, nCores = 1) %>% select(-row_derived_from_non_unique_identifiers)
  
  expect_equal(new, out2)
})

test_that("all_ok3", code={
  datIn <- dat3
  
  new <- tpp2dCurveFit(data = datIn, nCores = 1) %>% select(-row_derived_from_non_unique_identifiers)
  
  expect_equal(new %>% mutate(pEC50_quality_check = as.numeric(pEC50_quality_check)), 
               out3 %>% mutate(pEC50_quality_check = as.numeric(pEC50_quality_check)), 
               tolerance = 1e-5)
})

test_that("idCol_updated", code={
  datIn <- dat2
  
  new <- tpp2dCurveFit(data = datIn, nCores = 1) %>% select(-row_derived_from_non_unique_identifiers)
  
  expect_equal(attr(new, "importSettings")$uniqueIdCol, "Protein_ID")
})

test_that(desc="no_uniqueIdCol", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$uniqueIdCol <- NULL
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})

test_that(desc="uniqueIdCol_not_character", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$uniqueIdCol <- numeric()
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})


test_that(desc="uniqueIdCol_not_found", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$uniqueIdCol <- "nonsense"
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})

test_that(desc="no_nonZeroCols", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$nonZeroCols <- NULL
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})

test_that(desc="nonZeroCols_not_character", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$nonZeroCols <- numeric()
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})


test_that(desc="nonZeroCols_not_found", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$nonZeroCols <- "nonsense"
  
  expect_error(
    tpp2dCurveFit(data = datIn)
  )
})


test_that(desc="warning_deprecated_fct_arg1", code={
  datIn <- dat2
  
  expect_warning(tpp2dCurveFit(data = datIn, configFile = NA))
  
})

test_that(desc="warning_deprecated_fct_arg2", code={
  datIn <- dat2
  
  expect_warning(tpp2dCurveFit(data = datIn, naStrs = NA))
  
})

test_that(desc="warning_deprecated_fct_arg3", code={
  datIn <- dat2
  
  expect_warning(tpp2dCurveFit(data = datIn, fcStr = NA))
  
})

test_that(desc="warning_deprecated_fct_arg4", code={
  datIn <- dat2
  
  expect_warning(tpp2dCurveFit(data = datIn, idVar = NA))
  
})

test_that(desc="warning_deprecated_fct_arg5", code={
  datIn <- dat2
  
  expect_warning(tpp2dCurveFit(data = datIn, nonZeroCols = NA))
  
})

test_that(desc="data_missing", code={
  
  expect_error(
    tpp2dCurveFit()
  )
})

