datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "randomESet.Rdata"))


test_that("pec50QualCheck", {
  lbd <- 0
  ubd <- 10
  pec50 <- c(-5, 0, 5, 10, 15, -Inf, Inf, NA)
  expected <- c("< xmin", 0, 5, 10, NA, "< xmin", NA, NA)
  returned <- TPP:::curveFitFctCCR_pEC50qualCheckCol(x=pec50, xmin=lbd, xmax=ubd)
  expect_equal(expected ,returned)
})