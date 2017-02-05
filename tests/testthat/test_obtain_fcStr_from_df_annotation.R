dataPath <- system.file("test_data", package="TPP")

# Load function input:
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_normResults2.rds")) # example input from an older experiment (12 rows)

# Start tests:
test_that("find_norm_fc_col", code={
  # if both options are given, pick the column prefix for normalized fold changes:
  new <- TPP:::obtain_fcStr_from_df_annotation(dat2)
  expect_equal(new, "norm_rel_fc_protein_")
})

test_that("find_unmod_fc_col", code={
  # if slot for normalized fcstr is NULL, pick column prefix for unmodified fold changes:
  datIn <- dat2
  attr(datIn, "importSettings")$fcStrNorm <- NULL # set field content to null
  new <- TPP:::obtain_fcStr_from_df_annotation(datIn)
  expect_equal(new, "rel_fc_protein_")
})

test_that("missing_attr_entry", code={
  # if slot for normalized fcstr is NULL, pick column prefix for unmodified fold changes:
  datIn <- dat2
  attr(datIn, "importSettings")["fcStrNorm"] <- NULL # remove field completely
  new <- TPP:::obtain_fcStr_from_df_annotation(datIn)
  expect_equal(new, "rel_fc_protein_")
})

test_that("missing_unmod_str_but_valid_norm_str", code={
  # if slot for normalized fcstr is NULL, pick column prefix for unmodified fold changes:
  datIn <- dat2
  attr(datIn, "importSettings")["fcStr"] <- NULL # remove field completely
  new <- TPP:::obtain_fcStr_from_df_annotation(datIn)
  expect_equal(new, "norm_rel_fc_protein_")
})

test_that("invalid_norm_FC_colname", code={
  # specified fcStr prefix is not correct because no such column exists in table:
  datIn <- dat2
  attr(datIn, "importSettings")["fcStrNorm"] <- "nonsense" # this string does not exist in column names
  expect_error(TPP:::obtain_fcStr_from_df_annotation(datIn))
})

test_that("invalid_unmod_FC_colname", code={
  # specified fcStr prefix is not correct because no such column exists in table:
  datIn <- dat2
  attr(datIn, "importSettings")$fcStr <- "nonsense" # this string does not exist in column names
  attr(datIn, "importSettings")["fcStrNorm"] <- NULL # remove field completely
  expect_error(TPP:::obtain_fcStr_from_df_annotation(datIn))
})

test_that("invalid_unmod_FC_colname_but_valid_norm_FC_colname", code={
  # specified fcStr prefix is not correct because no such column exists in table:
  datIn <- dat2
  attr(datIn, "importSettings")["fcStr"] <- NULL # remove field completely
  new <- TPP:::obtain_fcStr_from_df_annotation(datIn)
  expect_equal(new, "norm_rel_fc_protein_")
})

test_that("missing_setting_list", code={
  # specified fcStr prefix is not correct because no such column exists in table:
  datIn <- dat2
  attr(datIn, "importSettings") <- NULL # remove field completely
  expect_error(TPP:::obtain_fcStr_from_df_annotation(datIn))
})
