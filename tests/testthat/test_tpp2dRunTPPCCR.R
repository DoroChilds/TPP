load(system.file("example_data/2D_example_data/referenceCCRConfig.RData", package="TPP"))
load(system.file("example_data/2D_example_data/exampleRunCCRInput.RData", package="TPP"))

test_that("evalRunCCR", code={
  CCRresults <- tpp2dRunTPPCCR(configFile=exampleCCRConfig, dataTable=exampleRunCCRInput, idVar="unique_ID")
  expect_equal(CCRresults$passed_filter[1], TRUE)
  expect_equal(CCRresults$compound_effect[1], "destabilized")
})

test_that("evalRunCCRErr1", code={
  expect_error(tpp2dRunTPPCCR(configFile=exampleCCRConfig, dataTable=exampleRunCCRInput, idVar="noneSense"))
})

test_that("evalRunCCRErr2", code={
  expect_error(tpp2dRunTPPCCR(configFile=exampleCCRConfig, dataTable=exampleRunCCRInput, idVar="representative"))
})

test_that("evalRunCCRErr3", code={
  expect_error(tpp2dRunTPPCCR(dataTable=exampleRunCCRInput, idVar="representative"))
})

test_that("evalRunCCRErr4", code={
  expect_error(tpp2dRunTPPCCR(configFile=exampleCCRConfig, idVar="representative"))
})