data(panobinostat_2DTPP_smallExample)
filePath <- system.file("example_data", package="TPP")

test_that(desc="evalConfig", code={
  configTable <- tpp2dEvalConfigTable(panobinostat_2DTPP_config)
  expect_identical(configTable, panobinostat_2DTPP_config)
})

# test_that(desc="evalConfigTxt", code={
#   cfgPath <- file.path(filePath, "2D_example_data/panobinostat_ex_confg.txt")
#    configTable <- tpp2dEvalConfigTable(configTable = cfgPath) %>% mutate(Path = NULL)
#    expect_identical(configTable, panobinostat_2DTPP_config)
#  })

# test_that(desc="evalConfigCsv", code={
#   #   cfgPath <- file.path(filePath, "2D_example_data/panobinostat_ex_confg.csv")
#   configTable <- tpp2dEvalConfigTable(configTable=file.path(cfgPath))
#   expect_identical(configTable, panobinostat_2DTPP_config)
# })

test_that(desc="evalErConfig", code={
  expect_error(tpp2dEvalConfigTable(configTable=NULL))
})