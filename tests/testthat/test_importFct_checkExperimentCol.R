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
## function 'importFct_checkExperimentCol':
## ------------------------------------------------------------------------- ##
test_that(desc="'confgExpCol_allok", code={
  oldCol <- c("abc", "def", "123")
  newCol <- TPP:::importFct_checkExperimentCol(expCol=oldCol)
  expect_equal(newCol, oldCol)
})

test_that(desc="'confgExpCol_convertAlnum", code={
  oldCol <- c("ok1", "abc!", "def'", "123$", "ok2")
  newCol <- TPP:::importFct_checkExperimentCol(expCol=oldCol)
  expect_equal(newCol, gsub("([^[:alnum:]])", "_", oldCol))
})

test_that(desc="'confgExpCol_NULLcol", code={
  oldCol <- NULL
  expect_error(TPP:::importFct_checkExperimentCol(expCol=oldCol))
})