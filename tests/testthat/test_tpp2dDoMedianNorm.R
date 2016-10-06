data(panobinostat_2DTPP_smallExample)
load(system.file("example_data/2D_example_data/referenceFcData.RData", package="TPP"))
load(system.file("example_data/2D_example_data/referenceNormData.RData", package="TPP"))

test_that(desc="evalMedianNormHead", code={
  NormData2d <- tpp2dDoMedianNorm(configTable = panobinostat_2DTPP_config, 
                                  dataTable = headData2dFc,
                                  fcStr = "rel_fc_protein_") 
  expect_identical(NormData2d, headNormData2d)
})

test_that(desc="evalMedianNormTail", code={
  NormData2d <- tpp2dDoMedianNorm(configTable = panobinostat_2DTPP_config, 
                                  dataTable = tailData2dFc,
                                  fcStr = "rel_fc_protein_") 
  expect_identical(NormData2d, tailNormData2d)
})

test_that(desc="evalMedianNormErr", code={
  expect_error(tpp2dDoMedianNorm(configTable = panobinostat_2DTPP_config, 
                                 dataTable = headData2dFc,
                                  fcStr="noneSense")) 
})