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
## function 'importFct_checkConditions':
## ------------------------------------------------------------------------- ##
## test 1: correct condition input
test_that(desc='checkCondCorrect', code={
  
  condIn <- refConds
  lenIn <- NULL
  
  condOut <- TPP:::importFct_checkConditions(condInfo = condIn, expectedLength = lenIn)
  
  expect_equal(condIn, condOut)
})

## test 1.2: conditions are lower case 
test_that(desc='checkCondCorrect', code={
  
  condIn <- c("vehicle", "vehicle", "treatment", "treatment")
  lenIn <- NULL
  
  condOut <- TPP:::importFct_checkConditions(condInfo = condIn, expectedLength = lenIn)
  
  expect_equal(refConds, condOut)
})

## test 2: NULL input because column was not present in config table -> create NAs, produce message
test_that(desc='checkCondNull', code={
  
  condIn <- NULL
  lenIn <- 4
  
  condOut <- TPP:::importFct_checkConditions(condInfo = condIn, expectedLength = lenIn)
  
  expect_equal(condOut, rep(NA_character_, 4))
  
})

## test 3: Conditions named differently than 'Vehicle' and 'Treatment' -> stop
test_that(desc='checkCond_wrongChars', code={
  
  condIn <- c("A","B", "C", "D")
  lenIn <- 4
  
  expect_error(
    TPP:::importFct_checkConditions(condInfo = condIn, expectedLength = lenIn)
  )
})

## test 4: Condition column contains numbers instead of characters -> stop
test_that(desc='checkCond_noChars', code={
  
  condIn <- 1:4
  lenIn <- 4
  
  expect_error(
    TPP:::importFct_checkConditions(condInfo = condIn, expectedLength = lenIn)
  )  
})