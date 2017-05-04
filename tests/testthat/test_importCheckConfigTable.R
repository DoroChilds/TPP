## ------------------------------------------------------------------------- ##
## test for 1D data:
## ------------------------------------------------------------------------- ##
data(hdacTR_smallExample)
data(hdacCCR_smallExample)
refLabels <- c('126', '127L', '127H', '128L', '128H','129L', '129H','130L', '130H', '131L')
refTemps  <- c(67, 63, 59, 56, 53, 50, 47, 44, 41, 37)
refTempMat <- matrix(rep(refTemps,4), nrow=4, byrow=TRUE, dimnames = list(NULL, refLabels))
refNames <- c("Vehicle_1", "Vehicle_2", "Panobinostat_1", "Panobinostat_2")
refConds <- c("Vehicle", "Vehicle", "Treatment", "Treatment")
refComps <- c(ComparisonVT1="Panobinostat_1_vs_Vehicle_1",
              ComparisonVT2="Panobinostat_2_vs_Vehicle_2")

test_that(desc="confgCheck", code={
  cfgIn <- hdacTR_config
  typeIn <- "TR"
  
  confgList <- TPP:::importCheckConfigTable(infoTable = cfgIn, type = typeIn)
  refList   <- list(expNames = refNames,
                    expCond  = refConds,
                    files    = NULL,
                    compStrs = refComps,
                    labels   = refLabels,
                    tempMatrix = refTempMat)
  expect_equal(confgList, refList)
})

test_that(desc="confgCheck_tibble", code={
  cfgIn <- hdacTR_config %>% as.tbl()
  typeIn <- "TR"
  
  confgList <- TPP:::importCheckConfigTable(infoTable = cfgIn, type = typeIn)
  refList   <- list(expNames = refNames,
                    expCond  = refConds,
                    files    = NULL,
                    compStrs = refComps,
                    labels   = refLabels,
                    tempMatrix = refTempMat)
  expect_equal(confgList, refList)
})

test_that(desc="confgCheck_expColNULL", code={
  hdacTR_config$Experiment <- NULL
  expect_error(TPP:::importCheckConfigTable(infoTable=hdacTR_config, type="TR"))
  expect_error(TPP:::importCheckConfigTable(infoTable=hdacTR_config, type="2D"))
})

test_that(desc="'confgCheck_expColConvertAlnum", code={
  oldCol <- paste(hdacTR_config$Experiment, c("", ":", "'", "!"))
  hdacTR_config$Experiment <- oldCol
  confgList <- TPP:::importCheckConfigTable(infoTable=hdacTR_config, type="TR")
  newCol <- confgList$expNames
  expect_equal(newCol, gsub("([^[:alnum:]])", "_", oldCol))
})

test_that(desc="'confgCheck_expColEmptyEntries", code={
  hdacTR_config$Experiment <- gsub("Vehicle_1", "", hdacTR_config$Experiment)
  confgList <- TPP:::importCheckConfigTable(infoTable=hdacTR_config, type="TR")
  expect_equal(confgList$expNames, hdacTR_config$Experiment[-1])
})

## ------------------------------------------------------------------------- ##
## test for 2D data:
## ------------------------------------------------------------------------- ##

filePath <- system.file("test_data", package="TPP")
load(file.path(filePath, "panobinostat_2DTPP_smallExample.RData"))

cfg <- panobinostat_2DTPP_config  %>% select(-Path)
cfgRef <- panobinostat_2DTPP_config %>% select(-Path)

test_that(desc="allOK", code={
  cfgIn <- cfg
  out <- TPP:::importCheckConfigTable(infoTable = cfgIn, type = "2D")
  expect_identical(cfgRef, out)
})

test_that(desc="evalConfigSpecialChars", code={
  cfgIn <- cfg %>% mutate_("126" = "'_'")
  ref <- cfgIn %>% mutate_("126" = "NULL")
  new <- TPP:::importCheckConfigTable(infoTable = cfgIn, type = "2D")
  expect_identical(ref, new)
})

test_that(desc="evalConfigTxt", code={
  cfgPath <- file.path(filePath, "panobinostat_ex_confg.txt")
  ref <- cfgRef
  new <- TPP:::importCheckConfigTable(infoTable = cfgPath, type = "2D")
  expect_identical(ref, new)
})


test_that(desc="configTable_NULL", code={
  expect_error( 
    TPP:::importCheckConfigTable(infoTable = NULL)
    )
})

test_that(desc="configTable_missing", code={
  expect_error( 
    TPP:::importCheckConfigTable()
    )
})


test_that(desc="configTable_1column", code={
  cfgIn <- cfg %>% select(Experiment)
  expect_error( 
    TPP:::importCheckConfigTable(infoTable = cfgIn, type = "2D")
    )
})

test_that(desc="configTable_noExperimentcolumn", code={
  cfgIn <- cfg %>% select(-Experiment)
  expect_error( 
    TPP:::importCheckConfigTable(infoTable = cfgIn, type = "2D")
    )
})

test_that(desc="invalid_file_name", code={
  cfgIn <- getwd()
  expect_error(
    TPP:::importCheckConfigTable(infoTable = cfgIn, type = "2D")
  )
})

test_that(desc="invalid_type", code={
  cfgIn <- cfg
  expect_error(
    TPP:::importCheckConfigTable(infoTable = cfgIn, type = "dummy")
  )
})

