qcPlotFctVennTable <- function(plotIDs, allIDs){
  tableDF = data.frame()
  for(name in names(allIDs)){
    lmntsIn <- length(plotIDs[[name]])
    lmntsTotal <- length(allIDs[[name]])
    lmntsOut <- lmntsTotal - lmntsIn
    row_df = data.frame(condition=name, passed_filter=lmntsIn, 
                        filtered_out=lmntsOut, Total=lmntsTotal)
    tableDF = rbind(tableDF, row_df)
  }
  return(tableDF)
}  

