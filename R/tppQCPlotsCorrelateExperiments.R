#' @title Plot pairwise relationships between the proteins in different TPP
#'   experiments.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' # Quality control (QC) plots BEFORE normalization:
#' tppQCPlotsCorrelateExperiments(tppData=tpptrData, 
#' annotStr="Non-normalized Fold Changes", ggplotTheme=ggplotTheme)
#' # Quality control (QC) plots AFTER normalization:
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=tpptrDefaultNormReqs())
#' tpptrDataNormalized <- tpptrNorm$normData
#' tppQCPlotsCorrelateExperiments(tppData=tpptrDataNormalized, 
#' annotStr="Normalized Fold Changes", ggplotTheme=ggplotTheme)
#' 
#' @return List of plots for each experiment.
#' @param tppData List of expressionSets with data to be plotted.
#' @param annotStr String with additional information to be added to the plot.
#' @param path Location where to store resulting plot.
#' @param ggplotTheme ggplot theme for the created plots.
#' @seealso \code{\link{tppDefaultTheme}}
#' @export
tppQCPlotsCorrelateExperiments <- function(tppData, annotStr="", path=NULL,
                                           ggplotTheme=tppDefaultTheme()){
  if (length(tppData) > 1){
    expNames <- names(tppData)
    plotObjList <- list()
    for (n1 in 1:(length(expNames)-1)){
      for (n2 in (n1+1):length(expNames)){
        expName1 <- expNames[n1]
        expName2 <- expNames[n2]
        dat1 = tppData[[expName1]]
        dat2 = tppData[[expName2]]
        
        ## Remove first temperature entry (only contains ones anyway)
        dat1 <- dat1[,dat1$temperature>min(dat1$temperature, na.rm=TRUE)]
        dat2 <- dat2[,dat2$temperature>min(dat2$temperature, na.rm=TRUE)]
        
        ## Only regard proteins that were detected in both experiments:
        commonProteins <- intersect(featureNames(dat1), featureNames(dat2))
        commonDat1 <- dat1[commonProteins,]
        commonDat2 <- dat2[commonProteins,]
        fc1 <- data.frame(exprs(commonDat1))
        fc2 <- data.frame(exprs(commonDat2))
        df1 <- reshape(data=fc1, idvar="protID", ids=row.names(fc1), 
                       times=paste(commonDat1$temperature, "\U00B0 C"), 
                       timevar="temperature",
                       varying=list(colnames(fc1)), direction="long", 
                       v.names="fc1")
        df2 <- reshape(data=fc2, idvar="protID", ids=row.names(fc2), 
                       times=paste(commonDat2$temperature, "\U00B0 C"), 
                       timevar="temperature",
                       varying=list(colnames(fc2)), direction="long", 
                       v.names="fc2")
        dfPlot <- merge(df1, df2)
        xy = geom_abline(intercept=0, slope=1, linetype=1, colour="red")
        p <- ggplot(dfPlot, aes_string(x="fc1", y="fc2"))
        p <- p + geom_point(alpha=0.2, na.rm=TRUE) + xy
        p <- p + scale_x_continuous(limits = c(0, 1.5)) + 
          scale_y_continuous(limits = c(0, 1.5))
        p <- p + xlab(expName1) + ylab(expName2) + ggtitle(paste(annotStr))
        p <- p + facet_wrap("temperature")
        
        plotName <- paste("QC_plots", expName1, "vs", expName2, sep="_")
        plotObjList[[plotName]] <- p
        ## Save plot:
        if (!is.null(path)){
          if (!file.exists(path)) dir.create(path, recursive=TRUE)      
          plotFile <- paste(plotName, "_", annotStr, ".pdf", sep="")
          ggsave(filename=file.path(path, plotFile), plot=p, width=20, 
                 height=25, units="cm")
        }
      }
    }
    return(plotObjList)
  } else {
    warning("QC plots of experiment-wise correlations can only be created when more than one experiments are given.")
  }
  
}