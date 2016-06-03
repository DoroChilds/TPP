## ------------------------------------------------------------------------- ##
## function 'applyCoeffs':
## ------------------------------------------------------------------------- ##
test_that(desc="testApplyCoeffs", code={
  x <- matrix(rnorm(10), nrow=2, ncol=5)
  eSet <- ExpressionSet(assayData=x)
  
  normCoeffs <- rnorm(5)
  eSetNormed <- applyCoeffs(eSet, normCoeffs)
  
  xNew <- unname(exprs(eSetNormed))
  xRef <- rbind(x[1,] * normCoeffs, x[2,] * normCoeffs)
  expect_equal(xNew, xRef)
})

test_that(desc = "testNormReqs1", code = {
  data("hdacTR_smallExample")
  tpptrData <- suppressMessages(tpptrImport(hdacTR_config, hdacTR_data))
  
  # 1.: Can we leave the field 'otherRequirements' empty?
  reqs <- tpptrDefaultNormReqs()
  reqs$otherRequirements$thresholdLower <- -Inf
  ref <- suppressMessages(tpptrNormalize(data=tpptrData, normReqs=reqs))
  reqs$otherRequirements <- reqs$otherRequirements[0,]
  newResult <- suppressMessages(tpptrNormalize(data=tpptrData, normReqs=reqs))
  expect_equal(exprs(ref[["normData"]][[1]]), 
               exprs(newResult[["normData"]][[1]]))
})

test_that(desc = "testNormReqs2", code = {
  data("hdacTR_smallExample")
  tpptrData <- suppressMessages(tpptrImport(hdacTR_config, hdacTR_data))
  
  # 2.: Desired filter columns in 'otherRequirements' don't exist in the data
  reqs <- tpptrDefaultNormReqs()
  ref <- suppressMessages(tpptrNormalize(data=tpptrData, normReqs=reqs))
  newRow <- data.frame(colName = "dummy", 
                       thresholdLower = -Inf, 
                       thresholdUpper = Inf)
  reqs$otherRequirements <- rbind(reqs$otherRequirements,
                                  newRow)
  newResult <- suppressMessages(tpptrNormalize(data=tpptrData, normReqs=reqs))
  expect_equal(exprs(ref[["normData"]][[1]]), 
               exprs(newResult[["normData"]][[1]]))
})
