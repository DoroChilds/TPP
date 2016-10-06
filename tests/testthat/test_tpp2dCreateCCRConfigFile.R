data(panobinostat_2DTPP_smallExample)
load(system.file("example_data/2D_example_data/referenceCCRConfig.RData", package="TPP"))

test_that("createCCRConfig", code={
  CCRconfig <- tpp2dCreateCCRConfigFile(configTable = panobinostat_2DTPP_config)
  expect_identical(CCRconfig, exampleCCRConfig)
})

test_that("createCCRConfigErr", code={
  expect_error(tpp2dCreateCCRConfigFile(configTable=NULL))
})