datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "randomESet.Rdata"))

test_that("CCR_normalization", {
  datNormalized <- tppccrNormalize(data=exampleSet)
  mediansNew <- apply(exprs(datNormalized), 2, median, na.rm=T)
  expect_equal(unname(mediansNew), rep(1, ncol(exampleSet)))
})