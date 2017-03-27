# Prepare function input:
data(hdacTR_smallExample)
data(hdacCCR_smallExample)
data(panobinostat_2DTPP_smallExample)

cfgCols_TR <- colnames(hdacTR_config)
cfgCols_CCR <- colnames(hdacCCR_config)
cfgCols_2D <- colnames(panobinostat_2DTPP_config)

allTMTcols <- c("126", "127L", "127H", "128L", "128H", "129L", "129H", "130L", 
                "130H", "131L")

test_that(desc="allOK_TR", code={
  
  colsIn <- cfgCols_TR
  
  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  
  expect_equal(colsOut, allTMTcols)
})

test_that(desc="allOK_CCR", code={
  
  colsIn <- cfgCols_CCR

  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  
  expect_equal(colsOut, allTMTcols)
})

test_that(desc="allOK_2D", code={
  
  colsIn <- cfgCols_2D
  
  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  
  expect_equal(colsOut, allTMTcols)
})

test_that(desc="allOK_replicateColumn", code={
  
  colsIn <- c(cfgCols_TR, "Replicate")
  
  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  
  expect_equal(colsOut, allTMTcols)
})

test_that(desc="allOK_numericLabels", code={
  
  colsIn <- c("Experiment", "0.01", "0.1")
  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  
  expect_equal(colsOut, c("0.01", "0.1"))
})

test_that(desc="colsMissing", code={
  
  expect_error(colsOut <- TPP:::detectLabelColumnsInConfigTable())
  
})

test_that(desc="colsNotCharacter", code={
  
  colsIn <- c(1:2, NA)
  colsOut <- TPP:::detectLabelColumnsInConfigTable(allColumns = colsIn)
  expect_equal(colsOut, colsIn)
  
})

