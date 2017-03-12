qcPlotFct_Bottleplot <- function(mpDiffs, minSlopes, isHit=NULL, strHit, strNoHit, 
                                 expName1, expName2, addHist, yLimVec){ 
  ## Generate QC plots to compare melting point differences to minimal slopes.

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
    histPlot <- histPlot + geom_histogram(data=plotDf, na.rm = TRUE,
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
    combinedPlot <- arrangeGrob(scatterPlot, histPlot, ncol=2, widths=c(3,1))
    titleObj <- textGrob(paste(expName1, "vs.", expName2), 
                         gp=gpar(fontsize=20,fontface='bold'), just="top")
    
    grid.arrange(combinedPlot, heights=c(1,0.1), main=titleObj)
  }
}

