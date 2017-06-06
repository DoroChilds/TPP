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
## function 'importCheckTemperatures':
## -------------------------------------------------------------------------- ##
test_that(desc = 'checkTemperatures_all_ok', code = {
  # Checking function checkTemperatures.
  ## Input:
  tmps <- hdacTR_config[,refLabels]
  ## Expected output:
  T_ref <- as.matrix(tmps)
  ## Test:
  T_new <- TPP:::importCheckTemperatures(temp = tmps)
  expect_equal(T_ref, T_new)  
})

test_that(desc = 'checkTemperatures_NA_values', code = {
  # Checking function checkTemperatures
  ## Input:
  tmps <- hdacTR_config[,refLabels]
  tmps[1:2,1:5] <- NA
  ## Expected output:
  T_ref <- as.matrix(tmps)
  ## Test:
  T_new <- TPP:::importCheckTemperatures(temp = tmps)
  expect_equal(T_ref, T_new)  
})

test_that(desc = 'checkTemperatures_NA_rows', code = {
  # Checking function checkTemperatures
  ## Input:
  tmps <- hdacTR_config[,refLabels]
  tmps[1:2,] <- NA
  ## Test:
  expect_error(TPP:::importCheckTemperatures(temp = tmps))
})