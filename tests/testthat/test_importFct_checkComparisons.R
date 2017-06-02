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
## function 'importFct_checkComparisons':
## ------------------------------------------------------------------------- ##
test_that(desc = 'checkComps_correct', code = {
  # Checking function importFct_checkComparisons.
  compStrs <- TPP:::importFct_checkComparisons(confgTable = hdacTR_config)
  refStrs <- c(ComparisonVT1 = "Panobinostat_1_vs_Vehicle_1",
               ComparisonVT2 = "Panobinostat_2_vs_Vehicle_2")
  expect_equal(compStrs, refStrs)
})

test_that(desc = 'checkComps_nullCol', code = {
  # Checking function importFct_checkComparisons.
  testtab <- hdacTR_config
  testtab$ComparisonVT1 <- testtab$ComparisonVT2 <- NULL
  compStrs <- TPP:::importFct_checkComparisons(confgTable = testtab)
  expect_equal(compStrs, NULL)
})