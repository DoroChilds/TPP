datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "randomESet.Rdata"))

test_that("CCR_normalization_allOK", {
  datNormalized <- tppccrNormalize(data=list("Experiment_1"=exampleSet))
  mediansNew <- apply(Biobase::exprs(datNormalized[[1]]), 2, median, na.rm=T)
  expect_equal(unname(mediansNew), rep(1, ncol(exampleSet)))
})

test_that("CCR_normalization_wrongInput", {
  expect_error(tppccrNormalize(data=exampleSet))
})


