dataPath <- system.file("test_data", package="TPP")

# Load function input:
dat1 <- readRDS(file.path(dataPath, "panobinostat_2D_fcResults.rds"))
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_fcResults2.rds")) # example input from an older experiment

# Load expected result for the given input:
out1 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults.rds"))
out2 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults2.rds"))

test_that(desc="all_ok1", code={
  datIn <- dat1
  
  new <- tpp2dNormalize(data = datIn)
  
  expect_equal(new, out1)
})

test_that(desc="all_ok2", code={
  datIn <- dat2
  
  new <- tpp2dNormalize(data = datIn)

  expect_equal(object = new, expected = out2)
})

test_that(desc="all_ok2_different_sorting", code={
  datIn <- dat2 %>% arrange(representative, temperature)
  
  new <- tpp2dNormalize(data = datIn)
  old <- out2 %>% arrange(representative, temperature)
  
  expect_equal(object = new, expected = old)
})

test_that(desc="all_ok2_scientificFotma", code={
  datIn <- dat2 %>% rename_("rel_fc_protein_1e-4" = "rel_fc_protein_1")
  
  new <- tpp2dNormalize(data = datIn)
  old <- out2 %>% rename_("rel_fc_protein_1e-4" = "rel_fc_protein_1",
                          "norm_rel_fc_protein_1e-4" = "norm_rel_fc_protein_1")
  
  expect_equal(object = new, expected = old)
})



test_that(desc="no_fcStr", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$fcStr <- NULL
  
  expect_error(
    tpp2dNormalize(data = datIn)
  )
  
})

test_that(desc="fcStr_not_character", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$fcStr <- numeric()
  
  expect_error(
    tpp2dNormalize(data = datIn)
  )
  
})


test_that(desc="fcStr_not_found", code={
  
  datIn <- dat2
  attr(datIn, "importSettings")$fcStr <- "nonsense"
  
  expect_error(
    tpp2dNormalize(data = datIn)
    )
  
})

test_that(desc="warning_deprecated_fct_arg1", code={
  datIn <- dat2
  
  expect_warning(tpp2dNormalize(data = datIn, configTable = NA))

})

test_that(desc="warning_deprecated_fct_arg2", code={
  datIn <- dat2
  
  expect_warning(tpp2dNormalize(data = datIn, fcStr = NA))
  
})

test_that(desc="data_missing", code={
  
  expect_error(
    tpp2dNormalize()
  )
})

