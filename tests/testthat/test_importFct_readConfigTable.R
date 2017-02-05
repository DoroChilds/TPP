dirExmpDat <- system.file("example_data", package="TPP")

data(hdacTR_smallExample)
data(hdacCCR_smallExample)
refLabels <- c('126', '127L', '127H', '128L', '128H','129L', '129H','130L', '130H', '131L')
refTemps  <- c(67, 63, 59, 56, 53, 50, 47, 44, 41, 37)
refTempMat <- matrix(rep(refTemps,4), nrow=4, byrow=TRUE, dimnames = list(c(1:4), refLabels))
refNames <- c("Vehicle_1", "Vehicle_2", "Panobinostat_1", "Panobinostat_2")
refConds <- c("Vehicle", "Vehicle", "Treatment", "Treatment")
refComps <- c(ComparisonVT1="Panobinostat_1_vs_Vehicle_1",
              ComparisonVT2="Panobinostat_2_vs_Vehicle_2")

## ------------------------------------------------------------------------- ##
## function 'importFct_readConfigTable':
## ------------------------------------------------------------------------- ##
test_that(desc='confgImportDF', code={
  # Checking function importFct_readConfigTable (import TR-config from data frame
  cfg <- TPP:::importFct_readConfigTable(cfg=hdacTR_config)
  expect_equal(cfg, hdacTR_config)
})

# test_that(desc='confgImportFileTR', code={
#   # Checking function importFct_readConfigTable (import TR-config from xlsx)
#   cfg <- importFct_readConfigTable(cfg=file.path(dirExmpDat, "TR_example_data", "Panobinostat_TPP-TR_config.xlsx"))
#   cfg$Path <- NULL
#   expect_equal(cfg, hdacTR_config)
# })

test_that(desc='confgImportFileCCR', code={
  # Checking function importFct_readConfigTable (import CCR-config from xlsx).
  cfg <- TPP:::importFct_readConfigTable(
    cfg=file.path(dirExmpDat, 
                  "CCR_example_data", 
                  "Panobinostat_TPP-CCR_config.xlsx")
    )
  cfg$Path <- NULL
  expect_equal(cfg, hdacCCR_config)
})