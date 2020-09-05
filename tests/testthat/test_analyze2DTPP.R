# Load function input:
data("panobinostat_2DTPP_smallExample")
cfg <- panobinostat_2DTPP_config
dat <- panobinostat_2DTPP_data

dataPath <- system.file("example_data", package="TPP") %>% file.path(., "2D_example_data")

# Start tests:
test_that(desc="allOK", code={
  
  cfgIn <- cfg
  datIn <- dat
  
  out <- analyze2DTPP(configTable = cfgIn,
                      data = datIn,
                      idVar = "representative",
                      addCol = c("clustername", "msexperiment_id"), 
                      intensityStr = "sumionarea_protein_",
                      qualColName = c("qupm", "qusm"),
                      nonZeroCols = "qusm",
                      fcStr = NULL)
  
  expect_equal(nrow(out), 4656)
})

test_that(desc="allOK_scientific_drug_concentration_format", code={
  ## Motivated by a bug report from Isabelle (2016-12-07)
  
  cfgIn <- cfg
  cfgIn[c(1,3,5,7,9,11), "128L"] <- "1E-4"
  cfgIn[c(1,3,5,7,9,11)+1, "130H"] <- "1E-4"
  
  datIn <- dat
  
  out<- analyze2DTPP(configTable = cfgIn,
                     data = datIn,
                     idVar = "representative",
                     addCol = c("clustername", "msexperiment_id"), 
                     intensityStr = "sumionarea_protein_",
                     qualColName = c("qupm", "qusm"),
                     nonZeroCols = "qusm",
                     fcStr = NULL)
  
  expect_equal(nrow(out), 4656)
  expect_equal(grep("1e.04", colnames(out), value = TRUE), 
                  c("norm_rel_fc_protein_1e.04_unmodified", 
                    "norm_rel_fc_protein_1e.04_normalized_to_lowest_conc",
                    "norm_rel_fc_protein_1e.04_transformed",
                    "sumionarea_protein_1e.04",
                    "rel_fc_protein_1e.04")) # The '-' character will be converted to '.' by the data.frame command when producing the wide result table.
  
})

test_that(desc="warning_deprecated_fct_arg", code={
  
  cfgIn <- cfg
  datIn <- dat
  
  expect_warning(
    analyze2DTPP(configFile = cfgIn,
                 data = datIn,
                 idVar = "representative",
                 addCol = c("clustername", "msexperiment_id"), 
                 intensityStr = "sumionarea_protein_",
                 qualColName = c("qupm", "qusm"),
                 nonZeroCols = "qusm",
                 fcStr = NULL)
  )
})

# test_that(desc="allOK_obtain_paths_from_config", code={
#   
#   cfgIn <- openxlsx::read.xlsx(file.path(dataPath, "Panobinostat_TPP-2D_config.xlsx")) %>% select(-Path)
#   
#   files <- tibble(Experiment = c("X020466", "X020467", "X020468", "X020469", "X020470", "X020471"),
#                   Path = c("PFC_39106_PF1_QV1_P0101868C.txt",
#                            "PFC_39108_PF1_QV1_P0101878C.txt",
#                            "PFC_39110_PF1_QV1_P0101888C.txt",
#                            "PFC_39101_PF1_QV1_P0101898C.txt",
#                            "PFC_39099_PF1_QV1_P0101908C.txt",
#                            "PFC_39093_PF1_QV1_P0101918C.txt") %>% file.path(dataPath, .))
#   
#   cfgIn <- cfgIn %>% left_join(files)
#   
#   out <- analyze2DTPP(configTable = cfgIn,
#                       idVar = "representative",
#                       addCol = c("clustername", "msexperiment_id"),
#                       intensityStr = "sumionarea_protein_",
#                       qualColName = c("qupm", "qusm"),
#                       nonZeroCols = "qusm",
#                       fcStr = NULL)
#   
#   check1 <- any(grepl("TPP_results", dir(dataPath))) # Was result path correctly created from the config file?
#   check2 <- all.equal(dim(out), c(61938, 46))
#   check3 <- all.equal(out$pEC50 %>% range(na.rm = TRUE), c(5.728183, 8.126123))
#   check4 <- all.equal(out$compound_effect %>% table %>% as.numeric, c(1921, 1092))
#   check5 <- all.equal(out$pEC50_outside_conc_range %>% table() %>% as.numeric(), c(1450, 1562))
#   expect_true(check1 & check2 & check3 & check4 & check5)
# })


## to do:
# test_that(desc="different_concentrations_per_experiment", code={ 
#   cfgIn <- cfg
#   cfgIn[1, "128L"] <- 0.01
#   
#   ## current to do: produce an informative error message stating that we currently can only handle the same concentrations across all temperatures.
#   ## long-term: enable different concentrations for different temperatures.
# })

