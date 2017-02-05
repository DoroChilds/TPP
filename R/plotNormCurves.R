plotNormCurves <- function(modelList, xMat, fcMat, r2Vec, nNormP, plotTheme){
  ## Plot QC plot to visualize curve fit of normalization curves to fold change
  ## medians.
  message("Creating QC plots to illustrate median curve fits.")
  grNames <- names(modelList)
  yMax <- 1.25
  theme_set(plotTheme)
  
  ## Prepare data objects for plotting:
  xLen      <- 100
  xMatLarge <- sapply(grNames, function(g) seq(from=min(xMat[,g]), 
                                               to=max(xMat[,g]), 
                                               length.out=xLen))
  yPred     <- sapply(grNames, function(g) {
    robustNlsPredict(model=modelList[[g]], newdata=list(x=xMatLarge[,g]))
    })
  
  plotDF_curves = data.frame()
  plotDF_points = data.frame()
  plotDF_anno <- data.frame()
  for(gn in grNames){
    plotDF_curves <- rbind(plotDF_curves, 
                          data.frame(Temperature = xMatLarge[, gn], 
                                     FoldChange = yPred[, gn],
                                     condition = rep(gn, nrow(xMatLarge))))
    
    plotDF_points <- rbind(plotDF_points, 
                          data.frame(Temperature = xMat[, gn], 
                                     FoldChange = fcMat[, gn],
                                     condition = rep(gn, nrow(xMat))))

    plotDF_anno <- rbind(plotDF_anno, 
                         data.frame(xPos = max(xMat)*0.95,
                                    yPos = yMax*0.95,
                                    lab = paste('R\U00B2', '=', 
                                                signif(r2Vec[gn], 3)),
                                    condition = gn ))
  }
  
  subtitle = paste("based on", nNormP, "proteins")
  
  p <- ggplot()
  p <- p + scale_y_continuous(limits = c(0, yMax))
  p <- p + theme(legend.position="none")
  p <- p + ggtitle(bquote(atop("Normalization curves", atop(.(subtitle)))))
  p <- p + xlab(paste('Temperature [\U00B0', 'C]', sep='')) + 
    ylab("Median fold change")
  p <- p + geom_line(data=plotDF_curves, na.rm = TRUE,
                     aes_string(x="Temperature", y="FoldChange", 
                                colour="condition"), size=1 )
  p <- p + geom_point(data=plotDF_points, na.rm = TRUE, size = 4,
                      aes_string(x="Temperature", y="FoldChange", 
                                 colour="condition"))
  p <- p + geom_text(data=plotDF_anno, 
                     aes_string(x="xPos", y="yPos", label="lab"))
  p <- p + facet_wrap(facets=~condition)
  
  return(p)
}