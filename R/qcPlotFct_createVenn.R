qcPlotFct_CreateVenn <- function(dataTab, grNames, compDF, mainTxt, method, 
                                 minR2){
  idList     <- list()
  filterList <- list()
  
  if (method == "individual"){    
    ## Create list of IDs per protein and apply filter to each of them:
    for (n in grNames){
      identifiedInExp <- dataTab[,paste("protein_identified_in", n, sep="_")]
      dataCurrentExp  <- dataTab[which(identifiedInExp),] 
      idList[[n]]     <- dataCurrentExp$Protein_ID
      
      ## filter all experiments according to R2:
      filterList[[n]] <- dataCurrentExp[,paste("R_sq", n, sep="_")] >= minR2      
    }
  } else if (method == "comparison"){
    for (i in 1:nrow(compDF)){
      ## Extract experiments for the current comparison
      expT  <- compDF$testGroup[i]
      expV  <- compDF$refGroup[i]
      compName <- paste(expT, expV, sep="_vs_")
      
      ## Extract protein IDs for the current comparison
      elemT <- dataTab[,paste("protein_identified_in", expT, sep="_")]
      elemV <- dataTab[,paste("protein_identified_in", expV, sep="_")]
      elemBoth <- elemT & elemV
      datBoth  <- dataTab[which(elemBoth),]
      idsBoth  <- datBoth$Protein_ID
      
      ## Retrieve the same filter that was used for p-value computation:
      filterColName <- paste("passed_filter", compName, sep="_")
      filterCol     <- datBoth[, filterColName]
      
      ## Store IDs and filters:
      idList[[compName]]     <- idsBoth
      filterList[[compName]] <- filterCol
      
    }
  }
  
  idListFiltered  <- sapply(names(idList), function(n) idList[[n]][which(filterList[[n]])], 
                            USE.NAMES=TRUE, simplify=FALSE)
  
  ## Produce plot:
  ## brewer.pal needs at least n=3 (bug?)
  plotCols <- brewer.pal(n = max(length(idListFiltered), 3), name = "Dark2")
  plotCols <- plotCols[1:length(idListFiltered)] ## If < 3 needed
  plotObj <- venn.diagram(x = idListFiltered, filename = NULL, fill = plotCols,
                          main = mainTxt, main.pos = c(0.5, 0.05),
                          main.cex = 1.5, main.fontfamily = "sans",
                          lty = 1, lwd = 0.25, alpha = 0.3, cex = 2,
                          label.col = 'blue', label.fontfamily = "sans",
                          cat.cex = 1.5, cat.pos = 0, cat.dist = (c(1:length(idListFiltered))/40)+0.05,
                          cat.fontface = 1, cat.fontfamily = "sans",
                          cat.col = plotCols)
  
  ## Add table that summarizes the proteins per group:  
  tableDF <- qcPlotFctVennTable(plotIDs=idListFiltered, allIDs=idList)
  if (method == "individual"){
    names(tableDF) <- gsub("passed_filter", paste("R2 >=", minR2), names(tableDF))
    names(tableDF) <- gsub("filtered_out" , paste("R2 <" , minR2), names(tableDF))
  } else if (method == "comparison") {
    # tableDF$condition <- gsub("_vs_", "_vs ", tableDF$condition)
  }
  plotObj <- addTableToPlot(plotObj = gTree(children=plotObj), meltVar = "condition",
                            tableDF = tableDF, clrs = plotCols)
  
  return(plotObj)
}

