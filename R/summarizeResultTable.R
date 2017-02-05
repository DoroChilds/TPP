summarizeResultTable <- function(inTable, wantedColPatterns, temperatures, 
                                 fcStr){
  
  #   lblsAndTemps = tppCfgTable[,c('Experiment', colnames(temperatures))]
  lbls = colnames(temperatures)
  colHeaders = names(inTable)
  
  ProteinIDColIdx = grep('Protein_ID', colHeaders)
  #clusterNameColIdx = grep('clustername', colHeaders)
  
  # annotBlocks = list()
  annotCols = data.frame(Protein_ID=inTable[, ProteinIDColIdx])
  # ctr = 0
  # for(cnci in clusterNameColIdx)
  # {
  #   ctr = ctr + 1
  #   tmp = data.frame('Protein_ID' = inTable[, ProteinIDColIdx],
  #                    'clustername' = inTable[, cnci])
  #   
  #   annotBlocks[[ctr]] = subset(tmp, !is.na(clustername))
  #   
  #   if(ctr == 1){
  #     annotCols = annotBlocks[[ctr]]
  #   } else {
  #     annotCols = merge(annotCols, annotBlocks[[ctr]], by=names(annotBlocks[[ctr]]), all=T)
  #   }
  # }
  
  outTableSummary = annotCols
  
  extraCols = c('qupm', 
                'ms1seqs', 
                'ms1intensity', 
                'ms1maxminusmin')
  
  for(ec in extraCols){
    extraColIdx <- grep(ec, colHeaders)
    #annotCols = cbind(annotCols, inTable[, extraColIdx])
    annotCols <- cbind(Protein_ID=inTable[, ProteinIDColIdx], inTable[, extraColIdx])
  }

  outTableDetail = annotCols
  
  cNames = names(outTableSummary)
  
  for(l in lbls){
    colIdx = grep(paste("norm_", fcStr, l, sep=""), colHeaders)
    
    tmp = inTable[,c(ProteinIDColIdx, colIdx)]
    outTableDetail = merge(outTableDetail, tmp, by="Protein_ID", all=TRUE)
    
    dataCols = c(2:ncol(tmp))
    tmpMedian = as.data.frame(cbind(tmp[,1], apply(tmp[,dataCols], 1, median, na.rm=TRUE)))
    names(tmpMedian) = c("Protein_ID", paste('median_FC', l, sep='_'))
    
    tmpSD = as.data.frame(cbind(tmp[,1], apply(tmp[,dataCols], 1, sd, na.rm=TRUE)))
    names(tmpSD) = c("Protein_ID", paste('sDev_FC', l, sep='_'))
  
    outTableSummary = merge(outTableSummary, tmpMedian, by="Protein_ID", all=TRUE)
    cNames = c(cNames, paste('median_FC', l, sep='_'))
    outTableSummary = merge(outTableSummary, tmpSD, by="Protein_ID", all=TRUE)
    cNames = c(cNames, paste('sDev_FC', l, sep='_'))
  }
  colnames(outTableSummary) <- cNames 
  
  cNames = names(outTableSummary)
  
  for(ptrn in wantedColPatterns){  
    colIdx = grep(ptrn, names(inTable))
    tmp = inTable[,c(ProteinIDColIdx, colIdx)]
    
    # outTableDetail = cbind(outTableDetail, tmp) 
    outTableDetail = merge(outTableDetail, tmp, by="Protein_ID", all=TRUE)
       
    dataCols = c(2:ncol(tmp))
    tmpMedian = as.data.frame(cbind(tmp[,1], apply(tmp[,dataCols], 1, median, na.rm=TRUE)))
    names(tmpMedian) = c("Protein_ID", paste('median', ptrn, sep='_'))
    
    tmpSD = as.data.frame(cbind(tmp[,1], apply(tmp[,dataCols], 1, sd, na.rm=TRUE)))
    names(tmpSD) = c("Protein_ID", paste('sDev', ptrn, sep='_'))
    
    outTableSummary = merge(outTableSummary, tmpMedian, by="Protein_ID", all=TRUE)
    cNames = c(cNames, paste('median', ptrn, sep='_'))
    
    outTableSummary = merge(outTableSummary, tmpSD, by="Protein_ID", all=TRUE)
    cNames = c(cNames, paste('sDev', ptrn, sep='_'))
    
  }
  colnames(outTableSummary) <- cNames
  
  return(list(detail=outTableDetail, summary=outTableSummary))
}
