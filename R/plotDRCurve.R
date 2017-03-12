plotDRCurve <- function(protID, fcDF, parDF, plotDir, allExp, addLegend, 
                        plotCols, verbose){
  ## Plot dose-response curves (one plot per protein).
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  experiment = param = concentration = foldChange = x = y <- NULL
  
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
  p <- p + geom_point(data=fcDF, na.rm = TRUE, size=4,
                      mapping=aes(x=concentration, y=foldChange, 
                                  colour=experiment))
  p <- p + geom_line(data=drCurves, na.rm = TRUE, mapping=aes(x=x, y=y, 
                                                              colour=experiment))
  p <- addTableToPlot(plotObj=p, tableDF=plotTable, meltVar = "condition", 
                      clrs=plotCols)  
  
  fName <- paste("drCurve_", gsub("([^[:alnum:]])", "_", protID),".pdf", 
                 sep="")
  fPath <- file.path(plotDir, fName)
  
  ## Print plot to PDF:
  pdf(file=fPath, width=7.87, height=9.84, useDingbats=FALSE)
  grid.arrange(p)
  dev.off()      
  ## We used to plot with the ggsave function until a major update in the 
  ## gridExtra package (to version 2.0.0) made ggsave incompatible with 
  ## gridArrange output. We'll later try to switch back to the command:
  # ggsave(filename=fPath, plot=p, width=20, height=25, units="cm", 
  # useDingbats=FALSE)
  
  return(data.frame("Protein_ID" = protID, "path"=fName))
}