# Load function input:
data("panobinostat_2DTPP_smallExample")
cfg <- panobinostat_2DTPP_config
dat <- panobinostat_2DTPP_data


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
  
  check1 <- nrow(out) == 4656
  
  expect_true(check1)
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
  
  check1 <- nrow(out) == 4656
  check2 <- all(grep("1e.04", colnames(out), value = TRUE) == 
                  c("norm_rel_fc_protein_1e.04_unmodified", 
                    "norm_rel_fc_protein_1e.04_normalized_to_lowest_conc",
                    "norm_rel_fc_protein_1e.04_transformed",
                    "sumionarea_protein_1e.04",
                    "rel_fc_protein_1e.04")) # The '-' character will be converted to '.' by the data.frame command when producing the wide result table.
  
  expect_true(check1 & check2)
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


## to do:
# test_that(desc="different_concentrations_per_experiment", code={ 
#   cfgIn <- cfg
#   cfgIn[1, "128L"] <- 0.01
#   
#   ## current to do: produce an informative error message stating that we currently can only handle the same concentrations across all temperatures.
#   ## long-term: enable different concentrations for different temperatures.
# })

