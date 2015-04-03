datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "validationData_meltCurveFilter.Rdata"))

test_that("missingFilterValues", {
  passedTest <- resultFilterCurvePars(r2V=refR2V, r2T=refR2T, plV=refPlV)
  newFilteredOut <- !passedTest
  expect_equal(sum(is.na(newFilteredOut)), sum(is.na(refFilteredOut)))
})

test_that("nonmissingFilterValues", {
  passedTest <- resultFilterCurvePars(r2V=refR2V, r2T=refR2T, plV=refPlV)
  newFilteredOut <- !passedTest
  naProteins <- is.na(newFilteredOut)
  expect_equal(refFilteredOut[!naProteins], newFilteredOut[!naProteins])
})
