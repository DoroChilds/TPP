# ---------------------------------------------------------#
# Unit tests for function 'pValFctPerformSingleComparison' #
# ---------------------------------------------------------#
datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "validationData_pvalComputation.Rdata"))

test_that("Invoke p-value computation for a single comparison (should work)", {
  i <- !refFilteredOut
  newP <- pValFctPerformSingleComparison(minsl=refMinSl[which(i)],
                                               mpdiff=refMpDiff[which(i)], 
                                               method="robustZ", 
                                               control=list(binWidth=300),
                                               comparisonName = "T_vs_V")
  expect_equal(newP, refP[i]*2, label = "Computed p-values")
})

test_that("No valid melting point difference", {
  i <- !refFilteredOut
  suppressWarnings(
    newP <- pValFctPerformSingleComparison(minsl=refMinSl[which(i)], 
                                           mpdiff=rep(NA, length(refMinSl)), 
                                           method="robustZ", 
                                           control=list(binWidth=300),
                                           comparisonName = "T_vs_V")
  )

  expect_equal(newP, rep(NA_real_, length(refMinSl)), label = "Computed p-values")
  
})

# test_that("Try different method (currently not supported)", {
#   i <- !refFilteredOut
#   newP <- try(pValFctPerformSingleComparison(minsl=refMinSl[which(i)], 
#                                              mpdiff=refMpDiff[which(i)], 
#                                              method="dummy", 
#                                              control=list(binWidth=300),
#                                              comparisonName = "T_vs_V"))
#   expect_error(newP)
# })

test_that("Given binWidth too big", {
  i <- !refFilteredOut
  expect_warning(
    pValFctPerformSingleComparison(minsl=refMinSl[which(i)], 
                                         mpdiff=refMpDiff[which(i)], 
                                         method="robustZ", 
                                         control=list(binWidth=3500),
                                         comparisonName = "T_vs_V")
  )
})
