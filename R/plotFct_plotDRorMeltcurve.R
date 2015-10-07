plotMeltingCurve <- function(modelList, listUpper=NULL, listLower=NULL,
                             xMat, fcMat, curvePars, protID, 
                             filename, plotTheme, expConditions, expComps, 
                             addLegend, useCI){
  
  # retrieving options via getOptions doesn't work with parallel execution
  # therefore options are passed as variables
  
  ## Produce a ggplot to vizualise fitted melting curves and their parameters.
  if(all(is.na(fcMat))){
    return(NULL)
    # return(FALSE)
  } else {
    grNames <- names(modelList)
    yMax <- 1.5
    theme_set(plotTheme)
    
    ## Prepare data objects for plotting:
    ## 1) Create long data tables for ggplot function:
    ## 1.1) Predicted values:
    xLen      <- 100
    xMatLarge <- sapply(grNames, function(g) 
      seq(from=min(xMat[g,]), to=max(xMat[g,]), length.out=xLen))
    yPred     <- sapply(grNames, function(g) 
      robustNlsPredict(model=modelList[[g]], newdata=list(x=xMatLarge[,g])))
    
    if (useCI) {
      
      yUpper <- sapply(grNames, function(g) 
        robustNlsPredict(model=listUpper[[g]], newdata=list(x=xMatLarge[,g])))
      
      yLower <- sapply(grNames, function(g) 
        robustNlsPredict(model=listLower[[g]], newdata=list(x=xMatLarge[,g])))
      
    }
    
    ## Determine order of group factors so that they form alternating
    ## Treatment/Vehicle pairs in the ggplot legend:
    compNumber <- assignCompNumber_to_expName(compDF=expComps, expNames=grNames)    
    plotCols   <- plotColors(expConditions, compNumber)
    plotlTypes <- plotLineTypes(expConditions)
    
    names(plotCols) <- grNames
    names(plotlTypes) <- grNames
    grOrder <- c()
    if (all(!is.na(compNumber))){
      for (r in unique(compNumber)){
        iT <- which(expConditions=="Treatment" & compNumber==r)
        iV <- which(expConditions=="Vehicle" & compNumber==r)
        grOrder <- c(grOrder, grNames[iT], grNames[iV])
      }      
    }
    if (length(grOrder) == 0) grOrder <- grNames
    plotCols <- plotCols[grOrder]
    plotlTypes <- plotlTypes[grOrder]
    
    groupCol1 <- factor(rep(grNames, each=xLen), levels=grOrder)
    plotData1 <- data.frame(Group       = groupCol1,
                            Temperature = numeric(length(grNames)*xLen),
                            FoldChange  = numeric(length(grNames)*xLen),
                            CiUp = rep(NA_real_, length(grNames)*xLen),
                            CiLow = rep(NA_real_, length(grNames)*xLen),
                            DataType    = "Model")
    for (g in grNames) {
      if (useCI) {
        plotData1[plotData1$Group==g, 
                  c("Temperature","FoldChange", "CiUp", "CiLow")] = cbind(xMatLarge[,g], yPred[,g], yUpper[,g], yLower[,g])
      } else {
        plotData1[plotData1$Group==g, 
                  c("Temperature","FoldChange")] = cbind(xMatLarge[,g], yPred[,g])
      }
    }
    
    ## 1.2) Measured values:
    groupCol2 <- factor(rep(grNames, each=ncol(fcMat)), levels=grOrder)
    plotData2 <- data.frame(Group       = groupCol2,
                            Temperature = numeric(length(grNames)*ncol(fcMat)),
                            FoldChange  = numeric(length(grNames)*ncol(fcMat)),
                            DataType    = "Measured")
    for (g in grNames) {
      plotData2[plotData2$Group==g, 
                c("Temperature","FoldChange")] = cbind(xMat[g,], fcMat[g,])
    }
    
    ## 1.3) Melting points:
    xMP <- subset(curvePars, select=meltPoint)
    yMP <- sapply(grNames, function(g) robustNlsPredict(modelList[[g]], 
                                                        newdata=list(x=xMP[g,])))
    
    groupCol3 <- factor(names(yMP), levels=grOrder)
    plotData3 <- data.frame(Group=groupCol3, yMP=yMP, xMP=xMP[,"meltPoint"])
    
    ## 2) Data frame with curve parameters:
    tableDF <- data.frame(condition = factor(grNames, levels=grOrder),
                          meltPoint = round(curvePars[grNames,"meltPoint"],2),
                          slope     = signif(curvePars[grNames,"slope"],2),
                          plateau   = round(curvePars[grNames,"plateau"],2),
                          R2        = round(curvePars[grNames,"R_sq"],2))
    
    ## Create plot:
    p <- ggplot()
    p <- p + scale_color_manual(values = plotCols)
    p <- p + scale_linetype_manual(values = plotlTypes)
    p <- p + scale_y_continuous(limits = c(0, yMax))
    if (addLegend){
      p <- p + theme(legend.position=c(1, 1), legend.justification=c(1,1), 
                     legend.title=element_blank())      
    } else {
      p <- p + theme(legend.position="none")      
    }
    p <- p + ggtitle(protID)
    p <- p + xlab(paste('Temperature [\U00B0', 'C]', sep='')) + 
      ylab("Fraction non-denatured")
    
    if (useCI) {
      p <- p + geom_ribbon(data=plotData1,
                           aes_string(x="Temperature", ymax="CiUp", ymin = "CiLow", 
                                      fill="Group"), alpha = 0.3) + 
        scale_fill_manual(values=plotCols)
    }
    
    p <- p + geom_line(data=plotData1, size=1,
                       aes_string(x="Temperature", y="FoldChange", 
                                  colour="Group", linetype="Group"))

    
    p <- p + geom_point(data=plotData2, na.rm=TRUE,
                        aes_string(x="Temperature", y="FoldChange", 
                                   colour="Group"))
    p <- p + geom_point(data=plotData3, shape=4, size=5, 
                        show_guide=FALSE, na.rm=TRUE,
                        aes_string(x="xMP", y="yMP", colour="Group"))
    
    p <- addTableToPlot(plotObj=p, tableDF=tableDF, meltVar="condition", 
                        clrs=plotCols)
   
    # print(p)
     
    ## Print plot to PDF:
    #if (!is.null(filename)) ggsave(filename=filename, plot=p, width=20, 
    #                               height=25, units="cm")    
    return(p)
  }
}

plotDRCurve <- function(protID, fcDF, parDF, plotDir, allExp, addLegend, plotCols, verbose){
  if (verbose) {
    message("Creating plot for protein ", protID)
  }
  yMax <- 1.5
  dose <- sort(unique(fcDF$concentration))
  
  ## Calculate DR curve values:
  xlim <- c(dose[2], dose[length(dose)])
  xvals <- seq(xlim[1],xlim[2], 0.01)
  drCurves <- c()
  plotTable <- c()
  for (e in allExp){
    if (any(parDF$experiment == e)){
      parsTmp <- subset(parDF, experiment==e)
      pec50 <- subset(parsTmp, param=="pEC50")$value
      slope <- subset(parsTmp, param=="slope")$value
      r2    <- subset(parsTmp, param=="R_sq")$value
      yPred <- eval(parse(text=fctSigmoidCCR()), 
                    envir=list(x=xvals, infl=-pec50, hill=slope))      
    } else {
      ## Add NA entries, if an experiment is missing for the current protein:
      pec50 = slope = r2 = NA
      yPred <- rep(NA_real_, length(xvals))
    }
    drCurves <- rbind(drCurves, data.frame(x=xvals, y=yPred, experiment=e))
    plotTable <- rbind(plotTable, data.frame(condition=as.factor(e), 
                                             pEC50 = signif(pec50,3),
                                             slope = signif(slope,3),
                                             R2    = signif(r2,3)))
  }
  
  ylim <- c(-0.5,yMax)
  p <- ggplot()
  p <- p + scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)
  p <- p + ggtitle(protID)
  p <- p + scale_color_manual(values = plotCols)
  if (addLegend){
    p <- p + theme(legend.position=c(1, 1), legend.justification=c(1,1), 
                   legend.title=element_blank())      
  } else {
    p <- p + theme(legend.position="none")      
  }
  p <- p + xlab("cpd. conc. (log M)") + ylab("normalized apparent stability" )
  p <- p + geom_point(data=fcDF, na.rm=TRUE, size=4,
                      mapping=aes(x=concentration, y=foldChange, colour=experiment))
  p <- p + geom_line(data=drCurves, mapping=aes(x=x, y=y, colour=experiment))
  p <- addTableToPlot(plotObj=p, tableDF=plotTable, meltVar = "condition", clrs=plotCols)  
  
  fName <- paste("drCurve_", gsub("([^[:alnum:]])", "_", protID),".pdf", sep="")
  fPath <- file.path(plotDir, fName)
  
  ## Print plot to PDF:
  pdf(file=fPath, width=7.87, height=9.84, useDingbats=FALSE)
  grid.arrange(p)
  dev.off()      
  ## We used to plot with the ggsave function until a major update in the 
  ## gridExtra package (to version 2.0.0) made ggsave incompatible with 
  ## gridArrange output. We'll later try to switch back to the command:
  # ggsave(filename=fPath, plot=p, width=20, height=25, units="cm", useDingbats=FALSE)
  
  return(data.frame("Protein_ID" = protID, "path"=fName))
}