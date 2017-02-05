dataPath <- system.file("test_data", package="TPP")

# Load function input:
dat1 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults.rds")) # example input
dat2 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults2.rds")) # example input from an older experiment (12 rows)
dat3 <- readRDS(file.path(dataPath, "panobinostat_2D_fitResults3.rds")) # example input from an older experiment (20 rows)

# Define output location:
outPath <- getwd()

# Start tests:
test_that("preserveCols", code={
  
  datIn <- head(dat1)
  
  f1 <- tpp2dExport(tab = datIn, outPath = outPath, 
                    addCol = NULL, trRef = NULL, addPlotColumns = FALSE)

  ref <- head(datIn[, !grepl("normalized_to_lowest_conc", colnames(datIn))])

  new <- openxlsx::read.xlsx(f1, sheet = "pEC50")
  file.remove(f1)
  
  refSorted <- ref[, colnames(new)]
  
  check1 <- all(sort(colnames(ref)) == sort(colnames(new)))
  check2 <- nrow(ref) == nrow(new)
  check3 <- all(sort(as.character(ref$representative)) == new$representative)
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="warning_deprecated_fct_arg1", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(configTable = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

test_that(desc="warning_deprecated_fct_arg2", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(resultPath = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

test_that(desc="warning_deprecated_fct_arg3", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(idVar = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

test_that(desc="warning_deprecated_fct_arg4", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(fcStr = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

test_that(desc="warning_deprecated_fct_arg5", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(intensityStr = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

test_that(desc="warning_deprecated_fct_arg6", code={
  datIn <- dat1
  
  expect_warning({
    f1 <- tpp2dExport(normalizedData = NA, tab = datIn, outPath = outPath)
    file.remove(f1)
    })
  
})

