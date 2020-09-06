dataPath <- system.file("test_data", package="TPP")

# Load function input:
datLeft <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults2.rds")) # example output from an older experiment (12 rows)

refPath <- file.path(system.file("data", package="TPP"),  "TPPTR_reference_results_HepG2.RData")

load(refPath)
datRight <- tppRefData$sumResTable$summary

# Start tests:
test_that("all_ok_reference_is_path", code={
  
  datIn <- datLeft
  refIn <- refPath
  
  new <- tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID")
  
  check1 <- all(sort(colnames(new)) == unique(sort(c(colnames(datLeft), colnames(datRight)))))
  check2 <- nrow(new) == nrow(datIn)
  check3 <- all.equal(datIn$representative, new$representative)
  
  expect_true(check1 & check2 & check3)
})

test_that("all_ok_reference_is_df", code={
  
  datIn <- datLeft
  refIn <- datRight
  
  new <- tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID")
  
  check1 <- all(sort(colnames(new)) == unique(sort(c(colnames(datLeft), colnames(datRight)))))
  check2 <- nrow(new) == nrow(datIn)
  check3 <- all.equal(datIn$representative, new$representative)
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="wrong_reference_path", code={
  
  datIn <- datLeft
  refIn <- "dummy"
  
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID")
  )
})


test_that(desc="idCols_mismatch_exists", code={
  
  datIn <- datLeft
  refIn <- datRight %>% filter(!(Protein_ID %in% c("IPI00013895.1", "IPI00030275.5")))
  
  expect_warning(new <- tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID"))

  missingRows <- new %>% 
    filter(representative %in% c("IPI00013895.1", "IPI00030275.5")) %>%
    subset(select = setdiff(colnames(datRight), "Protein_ID"))
  
  # The rows for the proteins that were missing in the reference data should be 
  # present, but filled with NAs:
  expect_equal(nrow(missingRows), 2)
  expect_true(all(is.na(missingRows)))
})


test_that(desc="idCols_no_match", code={
  
  datIn <- datLeft
  refIn <- datRight %>% filter(!(Protein_ID %in% datIn$representative))
  
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID")
  )
})


test_that(desc="IdCol_missing_2dData", code={
  
  datIn <- datLeft
  attr(datIn, "importSettings")$proteinIdCol <- "dummy"
  
  refIn <- datRight
  
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID")
  )
})

test_that(desc="IdCol_missing_refData", code={
  
  datIn <- datLeft
  refIn <- datRight
  
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "dummy")
  )
})

test_that(desc="warning_deprecated_fct_arg1", code={
  datIn <- datLeft
  refIn <- datRight
  
  expect_warning(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID",
                    idVar = NA)
  )
})

test_that(desc="warning_deprecated_fct_arg2", code={
  datIn <- datLeft
  refIn <- datRight
  
  expect_warning(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID",
                    data = NA)
  )
})


test_that(desc="warning_deprecated_fct_arg3", code={
  datIn <- datLeft
  refIn <- datRight
  
  expect_warning(
    tpp2dMerge2dRef(resultTable_2D = datIn, referenceDataSummary = refIn, refIDVar = "Protein_ID",
                    trRef = NA)
  )
})

test_that(desc="data_missing", code={
  expect_error(
    tpp2dMerge2dRef(referenceDataSummary = datRight, refIDVar = "Protein_ID")
  )
})

test_that(desc="refData_missing", code={
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datLeft, refIDVar = "Protein_ID")
  )
})


test_that(desc="refData_missing_invalid", code={
  expect_error(
    tpp2dMerge2dRef(resultTable_2D = datLeft, referenceDataSummary = NA, refIDVar = "Protein_ID")
  )
})

