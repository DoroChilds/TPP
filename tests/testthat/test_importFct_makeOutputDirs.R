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
## function 'importFct_makeOutputDirs':
## -------------------------------------------------------------------------- ##
test_that(desc = 'checkOutDirs_path_and_names_given', code = {
  ## Input:
  d <- getwd(); 
  f <- file.path(d, "dummy.txt")
  ## Expected output:
  e <- list(doWrite=TRUE, pathDataObj=file.path(d, "dataObj"), outDir=d)
  ## Test:
  res <- TPP:::importFct_makeOutputDirs(outDir=d, fNames=f)
  unlink(res$pathDataObj, recursive=TRUE)
  expect_equal(res, e)
})

test_that(desc = 'checkOutDirs_path_NULL', code = {
  fIn <- file.path(getwd(), "dummy.txt")
  
  res <- TPP:::importFct_makeOutputDirs(outDir=NULL, fNames = fIn)
  
  unlink(res$outDir, recursive=TRUE)
  expect_equal(dirname(res$outDir), dirname(fIn))
})

test_that(desc = 'checkOutDirs_files_NULL', code = {
  ## Input:
  d <- getwd()
  ## Expected output:
  e <- list(doWrite=TRUE, pathDataObj=file.path(d, "dataObj"), outDir=d)
  ## Test:
  res <- TPP:::importFct_makeOutputDirs(outDir=d, fNames=NULL)
  unlink(res$pathDataObj, recursive=TRUE)
  expect_equal(res, e)
})

test_that(desc = 'checkOutDirs_both_NULL', code = {
  ## Expected output:
  e <- list(doWrite=FALSE, pathDataObj=NULL, outDir=NULL)
  ## Test:
  res <- TPP:::importFct_makeOutputDirs(outDir=NULL, fNames=NULL)
  expect_equal(res, e)
})

