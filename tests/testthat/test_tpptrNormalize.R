data("hdacTR_smallExample")

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

reqs <- tpptrDefaultNormReqs()

test_that(desc = "testNormReqs1", code = {
  # 1.: Can we leave the field 'otherRequirements' empty?
  
  reqs$otherRequirements$thresholdLower <- -Inf
  
  ref <- suppressMessages(
    tpptrNormalize(data=tpptrData, normReqs=reqs)
    )
  
  reqs$otherRequirements <- reqs$otherRequirements[0,]
  
  newResult <- suppressMessages(
    tpptrNormalize(data=tpptrData, normReqs=reqs)
    )
  
  expect_equal(Biobase::exprs(ref[["normData"]][[1]]), 
               Biobase::exprs(newResult[["normData"]][[1]]))
})

test_that(desc = "testNormReqs2", code = {
  # 2.: Desired filter columns in 'otherRequirements' don't exist in the data
  
  ref <- suppressMessages(
    tpptrNormalize(data=tpptrData, normReqs=reqs)
    )
  
  newRow <- data.frame(colName = "dummy", 
                       thresholdLower = -Inf, 
                       thresholdUpper = Inf)
  
  reqs$otherRequirements <- rbind(reqs$otherRequirements, newRow)
  
  newResult <- suppressWarnings(
    tpptrNormalize(data=tpptrData, normReqs=reqs)
    )
  
  expect_equal(Biobase::exprs(ref[["normData"]][[1]]), 
               Biobase::exprs(newResult[["normData"]][[1]]))
})