datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "validationData_pvalComputation.Rdata"))

# test_that("proteinNumber",{
#   expect_equal(length(refP), 5732)
# })

test_that("minSlopes",{
  newMinSl <- computeMinimalSlopes(slV=refSlV, slT=refSlT)
  nonNa <- !is.na(newMinSl) & !is.na(refMinSl)
  expect_equal(newMinSl[nonNa], refMinSl[nonNa])
})

test_that("pvalue computation", {
  newMinSl <- computeMinimalSlopes(slV=refSlV, slT=refSlT)
  newMpDiff <- refMpT - refMpV
  i <- !refFilteredOut
  newP <- computePvalues(ids=refIDs[i], minSlopes=newMinSl[i], mpDiffs=newMpDiff[i], binWidth=300)
  expect_equal(newP, refP[i])
})

