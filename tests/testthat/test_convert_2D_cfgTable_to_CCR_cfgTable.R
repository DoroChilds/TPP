data("panobinostat_2DTPP_smallExample")

test_that("createCCRConfig", code={
  cfgIn <- panobinostat_2DTPP_config
  
  new <- TPP:::convert_2D_cfgTable_to_CCR_cfgTable(configTable = cfgIn)
  
  ref <- data.frame("Panobinostat", 5, 1, 0.143, 0.02, 0, row.names = "") %>%
    set_colnames(c("Experiment", 5, 1, 0.143, 0.02, 0))
  
  expect_equal(new, ref)
})

test_that("createCCRConfigErr", code={
  expect_error(
    TPP:::convert_2D_cfgTable_to_CCR_cfgTable(configTable=NULL)
    )
})
