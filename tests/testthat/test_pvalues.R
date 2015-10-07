datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "validationData_pvalComputation.Rdata"))

test_that("minSlopes",{
  newMinSl <- computeMinimalSlopes(xV=refSlV, xT=refSlT)
  nonNa <- !is.na(newMinSl) & !is.na(refMinSl)
  expect_equal(newMinSl[nonNa], refMinSl[nonNa])
})

test_that("mpDiffs",{
  newDiffs <- computeMPdiffs(xV=refMpV, xT=refMpT)
  expect_equal(newDiffs, refMpDiff)
})

test_that("binning", {
  message("Testing function 'assignBins'")
  load(file.path(datapath, "validationData_pvalBinning.Rdata"))
  newBins <- assignBins(ref_x, 10)
  expect_equal(newBins, refBins)
})

test_that("pvalue computation", {
  i <- !refFilteredOut
  newP <- computePvalues(minSlopes=refMinSl[i], mpDiffs=refMpDiff[i], binWidth=300)
  expect_equal(newP, refP[i])
})

