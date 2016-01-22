qcPlotFct_VennWrapper <- function(resultTab, expNames, grConditions, compDF, 
                                  minR2){  
  ## --------------------------------------------------------------------------
  ## I) Create a Venn diagram of protein numbers in all experiments
  ## --------------------------------------------------------------------------
  if(length(expNames) <= 5){
    message("Creating Venn diagram to summarize all experiments.")
    mainTxt = paste(nrow(resultTab), 'proteins identified overall')
    plotObj <- qcPlotFct_CreateVenn(dataTab=resultTab, grNames=expNames,
                                    mainTxt=mainTxt, method="individual", 
                                    minR2=minR2)
    grid.newpage()
    grid.draw(plotObj)
  } else {
    message("Venn diagrams can only be created for up to 5 experiments.")
  }
  ## --------------------------------------------------------------------------
  ## II) Creates separate Venn diagrams to compare the experiments in each 
  ## user-defined comparison pair (if applicable)
  ## --------------------------------------------------------------------------
  if (!is.null(compDF)){
    if(nrow(compDF) <= 5){
      message("Creating Venn diagram to summarize all comparisons.")
      plotObj <- qcPlotFct_CreateVenn(dataTab=resultTab, compDF=compDF, 
                                      mainTxt="", method="comparison")  
      grid.newpage()
      grid.draw(plotObj)
    } else {
      message("Venn diagrams can only be created for up to 5 comparisons.")
    }
  }
}

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

qcPlotFct_invokeBottleplots <- function(resultTable, compDF){
  alpha = 0.05 # significance level
  
  ## Retrieve experiment names and annotation:
  for(i in 1:nrow(compDF)){
    expNameV <- compDF[i, "refGroup"]
    expNameT <- compDF[i, "testGroup"]
    compName <- compDF[i, "name"]
    nameMpDiff    <- paste("diff_meltP", compName, sep="_")
    nameMinSlope  <- paste("min_slope",compName, sep="_")
    namePVal      <- paste("pVal_adj", compName, sep="_")
    
    xMpDiff <- resultTable[,nameMpDiff]
    xMinSl  <- resultTable[,nameMinSlope]
    xpVals  <- resultTable[,namePVal]
    isHit   <- xpVals <= alpha
    
    if (any(!is.na(isHit))){
      qcPlotFct_Bottleplot(mpDiffs=xMpDiff, minSlopes=xMinSl, isHit=isHit, 
                           strHit=paste("p_adj <=",alpha), 
                           strNoHit=paste("p_adj >",alpha),
                           expName1=expNameT, expName2=expNameV, addHist=TRUE, 
                           yLimVec=c(0, -1.5))
    } else {
      xStr1 <- paste('Melting point difference [\U00B0', 'C]', sep='')
      xStr2 <- paste('(',expNameT,' - ', expNameV,')', sep='')
      xLab <- paste(xStr1, '\n', xStr2, sep='')
      
      pBlank <- ggplot(data=data.frame(x=c(-1,-1,1,1), y=c(-1,1,-1,1)), aes(x=x,y=y)) + 
        geom_blank() + ylab("Minimal slope") + xlab(xLab) +
        geom_text(size=10, aes(x=0, y=0, label="No significant Tm shifts shown\nhere because p-value calculation\nwas not possible for any protein\nin this comparison.")) 
      print(pBlank)
    }
  }
  return(NULL)
}

qcPlotFct_Bottleplot <- function(mpDiffs, minSlopes, isHit=NULL, strHit, strNoHit, 
                                 expName1, expName2, addHist, yLimVec){ 
  ## Create 'bottle plot' and histogram of minimal slopes
  
  ## Create dataframe to be passed to ggplot function:
  plotDf <- data.frame(mpDiffs=mpDiffs, minSlopes=minSlopes)
  
  ## Prepare plot:
  xLabels <- c(-10, rep("",4), -5, rep("",4), 0, rep("",4), 5, rep("",4), 10)
  scatterPlot <- ggplot()
  scatterPlot <- scatterPlot + scale_alpha(range = c(0.3, 0.5))
  
  ## Only include proteins, for which p-values could be computed 
  ## (criterion: R2 > 0.8 (both columns) + Plateau < 0.3 (Vehicle)):
  if (is.null(isHit)){
    isHit <- rep(FALSE, nrow(plotDf))
    strHit <- "h"
    strNoHit <- "nh"
    addLegend <- FALSE
  } else {
    addLegend <- TRUE
  }
  
  naHits <- is.na(isHit)
  plotDf <- plotDf[!naHits, ]
  isHit <- isHit[!naHits]
  
  signfTest <- rep(NA, nrow(plotDf))
  signfTest[isHit]  <- strHit
  signfTest[!isHit] <- strNoHit
  plotDf$signfTest <- factor(signfTest, levels=c(strNoHit, strHit))
  scatterPlot <- scatterPlot + geom_point(data=plotDf, na.rm = TRUE,
                                          aes_string(x="mpDiffs", y="minSlopes", 
                                                     color="signfTest"), size=3)
  
  scatterPlot <- scatterPlot + scale_x_continuous(limits=c(-15,15),
                                                  breaks=seq(-10,10,by=1), 
                                                  labels=xLabels) 
  scatterPlot <- scatterPlot + ylim(yLimVec)
  xStr1 <- paste('Melting point difference [\U00B0', 'C]', sep='')
  xStr2 <- paste('(',expName1,' - ', expName2,')', sep='')
  xStr3 <- paste('n =', nrow(plotDf))
  xLab <- paste(xStr1, '\n', xStr2, '\n\n', xStr3, sep='')
  
  scatterPlot <- scatterPlot + xlab(xLab) + ylab("Minimal slope")
  if (addLegend){
    scatterPlot <- scatterPlot + theme(legend.position=c(0.85,0.85), 
                                       legend.title=element_blank())    
  } else {
    scatterPlot <- scatterPlot + theme(legend.position="none") 
  }
  scatterPlot <- scatterPlot + geom_point(data=subset(plotDf, signfTest==strHit), 
                                          na.rm = TRUE,
                                          aes_string(x="mpDiffs", y="minSlopes", 
                                                     color="signfTest"), size=3)
  if (!addHist){
    print(scatterPlot)
    return(scatterPlot)
  } else{
    histPlot <- ggplot()
    histPlot <- histPlot + geom_histogram(data=plotDf, 
                                          aes_string(x="minSlopes"), 
                                          binwidth=abs(diff(yLimVec))/100) 
    histPlot <- histPlot + coord_flip()
    if (yLimVec[1] < yLimVec[2]){
      histPlot <- histPlot + scale_x_continuous(labels=NULL, limits=yLimVec)            
    } else {
      histPlot <- histPlot + scale_x_reverse(labels=NULL, limits=yLimVec)      
    }
    histPlot <- histPlot + xlab(NULL) + ylab("Count\n\n\n")
    
    ## Save plot:
    combinedPlot <- suppressWarnings(arrangeGrob(scatterPlot, histPlot, ncol=2, 
                                                 widths=c(3,1)))
    titleObj <- textGrob(paste(expName1, "vs.", expName2), 
                         gp=gpar(fontsize=20,fontface='bold'), just="top")
    
    grid.arrange(combinedPlot, heights=c(1,0.1), main=titleObj)
  }
}

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

qcPlotFct_MeltPointHist <- function(resultTab=resultTable, expNames=expNames, 
                                    minR2, expConds){
  if (length(expNames) > 1){
    combis <- combn(expNames,2)
    numCompares <- ncol(combis)  
    
    diffDF = data.frame()
    for (j in 1:numCompares){
      n1 <- combis[1,j] 
      n2 <- combis[2,j]
      
      ## Compute melting point differences (ensure direction 'Treatment' -
      ## 'Vehicle', if information about conditions available)
      cond1 <- unname(expConds[expNames==n1])
      cond2 <- unname(expConds[expNames==n2])
      if (identical(cond2, "Treatment") && identical(cond1, "Vehicle")){
        n1_tmp <- n1
        n1 <- n2
        n2 <- n1_tmp
      }
      mpCol_minuend    <- paste("meltPoint", n1, sep="_")
      mpCol_subtrahend <- paste("meltPoint", n2, sep="_")
      
      mp1 <- resultTab[, mpCol_minuend]
      mp2 <- resultTab[, mpCol_subtrahend]
      mpDiff <- mp1 - mp2
      
      #  filter data according to R_sq
      r2_1 <- resultTab[, paste("R_sq", n1, sep="_")]
      r2_2 <- resultTab[, paste("R_sq", n2, sep="_")]
      mpDiff <- mpDiff[which(r2_1 >= 0.8 & r2_2 > 0.8)]
      
      # Compute median and standard deviation of differences
      mdn = round(median(mpDiff, na.rm=TRUE),3)
      q = quantile(mpDiff, probs = c(0.1587, 0.8413), na.rm=TRUE)
      robustSD_left  = round(q[1], 3)
      robustSD_right = round(q[2], 3)
      
      # Create plot annotation:
      titl=paste('median = ', mdn, ' , \n robust SD: left = ', robustSD_left, 
                 ' | right = ', robustSD_right, sep='')
      
      compName <- paste(n1, "-", n2)
      if (length(mpDiff)>0){ # do not produce plot in case no Tm diff exists for 
                             # the given combination (for example, because of 
                             # poor curve fits).

        
      tmpDF = data.frame('diff' = as.numeric(mpDiff), 'comparison'= compName)
      
      p = ggplot(data = tmpDF)
      p = p + stat_bin(aes(x=diff), colour='black', alpha=0.3, geom="bar", binwidth=0.2)
      xLabels <- seq(-10,10,by = 5)
      p = p + scale_x_continuous(limits=c(-15,15), breaks=xLabels, labels=xLabels)
      p = p + theme(legend.position="bottom")
      xStr1 <- paste('Melting point difference [\U00B0', 'C]', sep='')
      xStr2 <- paste('(',compName,')', sep='')
      xStr3 <- paste('n =', length(mpDiff), ', R_sq >=', minR2)
      p = p + xlab(paste(xStr1, '\n', xStr2, '\n\n', xStr3))
      p = p + ggtitle(titl)
      print(p)
      
      diffDF = rbind(diffDF, tmpDF)
      #       } else{ # to do!
      #         ggplot(data=data.frame(x=c(-1,-1,1,1), y=c(-1,1,-1,1)), aes(x=x,y=y)) + 
      #           geom_blank() + ylab("Minimal slope") + xlab(xLab) +
      #           geom_text(size=10, aes(x=0, y=0, label="No significant Tm shifts shown here\nbecause p-value calculation was not\npossible for any protein.")) 
      }
    }
    
    p = ggplot(data = diffDF)
    p = p + stat_bin(aes(x=diff, colour=comparison, fill=comparison), 
                     alpha = 0.05, geom="area", binwidth=0.2, 
                     position=position_dodge(width = 0))
    p = p + stat_bin(aes(x=diff, colour=comparison), 
                     geom="line", binwidth=0.2, position=position_dodge(width = 0))
    xLabels <- seq(-10,10,by = 5)
    p = p + scale_x_continuous(limits=c(-15,15), breaks=xLabels, labels=xLabels)
    p = p + theme(legend.position="bottom")
    p = p + guides(colour = guide_legend(nrow = max(1, floor(numCompares/2) ))) 
    p = p + xlab('melting point difference')
    print(p)
  } else {
    message("No melting point difference histograms were produced because they 
            require at least two experiments.")
  }
}
