extractConc <- function(configTable){
  # Original code from Nils:
  # exp.cols <- which(colnames(configTable) %in% (colnames(configTable)[-which(colnames(configTable) %in% 
  #                             c("Compound", "Experiment", "Temperature", "RefCol", "Path"))]))
  # conc.cols <- sapply(exp.cols, function(x) return(levels(as.factor(configTable[[x]]))[which(levels(as.factor(configTable[[4]]))!="-")]))
  # return(conc.cols)
  
  # # Same code, but with improved readability:
  # nonExpCols <- c("Compound", "Experiment", "Temperature", "RefCol", "Path")
  # allCols <- colnames(configTable)
  # labelCols <- setdiff(allCols, nonExpCols)
  # 
  # firstColName <- labelCols[1]
  # firstColLevels <- configTable[[firstColName]] %>% as.character %>% unique
  # firstColValid <- which(firstColLevels != "-")
  # 
  # labelColsPos <- match(labelCols, allCols)
  # conc.cols <- sapply(labelColsPos,
  #                     function(i){
  #                       colTmp <- configTable[[i]]
  #                       colTmpLevels <- unique(as.character(colTmp))
  #                       out <- colTmpLevels[firstColValid]
  #                       return(out)
  #                     })
  # return(conc.cols)
  
  # Own code:
  nonExpCols <- c("Compound", "Experiment", "Temperature", "RefCol", "Path")
  allCols <- colnames(configTable)
  labelCols <- setdiff(allCols, nonExpCols)

  labelValues <- configTable[,labelCols]
  labelValuesUnique <- labelValues %>% apply(2, unique)
  labelValuesNum <- labelValuesUnique[labelValuesUnique != "-"] %>% unname
  return(labelValuesNum)
}