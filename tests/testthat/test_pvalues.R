datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "validationData_pvalComputation.Rdata"))

test_that("minSlopes",{
  newMinSl <- TPP:::computeMinimalSlopes(xV=refSlV, xT=refSlT)
  nonNa <- !is.na(newMinSl) & !is.na(refMinSl)
  expect_equal(newMinSl[nonNa], refMinSl[nonNa])
})

test_that("mpDiffs",{
  newDiffs <- TPP:::computeMPdiffs(xV=refMpV, xT=refMpT)
  expect_equal(newDiffs, refMpDiff)
})

test_that("binning", {
  message("Testing function 'assignBins'")
  load(file.path(datapath, "validationData_pvalBinning.Rdata"))
  newBins <- TPP:::assignBins(ref_x, 10)
  expect_equal(newBins, refBins)
})

test_that("pvalue computation", {
  i <- !refFilteredOut
  newP <- TPP:::computePvalues(minSlopes = refMinSl[i], 
                         mpDiffs = refMpDiff[i], 
                         binWidth = 300,
                         type = "p", 
                         pAdj = "fdr")
  expect_equal(newP/2, refP[i])
})

