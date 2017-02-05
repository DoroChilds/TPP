testDataPath <- system.file("test_data", package="TPP")

# Load function input:
data("panobinostat_2DTPP_smallExample")
cfg <- panobinostat_2DTPP_config
dat <- panobinostat_2DTPP_data

# Load expected result for the given input:
results <- readRDS(file.path(testDataPath, "panobinostat_2D_importResults.rds"))

# Start tests:
test_that(desc="evalImportData", code={
  cfgIn <- cfg
  datIn <- dat
  ref <- results
  new <- tpp2dImport(configTable = cfgIn,
                     data = datIn,
                     idVar = "representative",
                     addCol = c("clustername", "msexperiment_id"), 
                     intensityStr = "sumionarea_protein_",
                     qualColName = c("qupm", "qusm"),
                     nonZeroCols = "qusm",
                     fcStr = NULL)
  # # The reference data were created before introduction of attributes in version 2.99.0:
  # attr(new, "importSettings")  <- NULL
  # attr(new, "configTable") <- NULL
  expect_identical(ref, new)
})

test_that(desc="cfg_empty_entries", code={
  ## After bug report from Isabelle (01.02.2017)
  cfgIn <- cfg %>% mutate(`126` = gsub("-", NA, `126`))
  datIn <- dat
  
  ref <- results
  attr(ref, "configTable") <- cfgIn %>% select(-Path)
  
  new <- tpp2dImport(configTable = cfgIn,
                     data = datIn,
                     idVar = "representative",
                     addCol = c("clustername", "msexperiment_id"), 
                     intensityStr = "sumionarea_protein_",
                     qualColName = c("qupm", "qusm"),
                     nonZeroCols = "qusm",
                     fcStr = NULL)
  # # The reference data were created before introduction of attributes in version 2.99.0:
  # attr(new, "importSettings")  <- NULL
  # attr(new, "configTable") <- NULL
  expect_identical(ref, new)
})



test_that(desc="evalImportDataErr1", code={
  configTable <- cfg
  data <- NULL
  expect_error(tpp2dImport(configTable = configTable, 
                               data = data, 
                               idVar = "representative",
                               intensityStr = "sumionarea_protein_",
                               addCol = c("qusm","clustername", "msexperiment_id"),
                               qualColName = "qupm", 
                               fcStr = NULL))
})

test_that(desc="evalImportDataErr2", code={
  configTable <- NULL
  data <- dat
  expect_error(tpp2dImport(configTable = configTable, 
                               data = data, 
                               idVar = "representative",
                               intensityStr = "sumionarea_protein_",
                               addCol = c("qusm","clustername", "msexperiment_id"),
                               qualColName = "qupm", 
                               fcStr = NULL))
})
