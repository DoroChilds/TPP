#' @title Venn diagrams of detected proteins per experiment.
#' @description \code{tppVenn} illustrates the overlaps between different 
#'   TPP-TR/CCR experiments by a venn diagram.
#' @param data list of ExpressionSets that contain the imported data per 
#'   experiment (return value of function \code{\link{tpptrImport}} or 
#'   \code{\link{tppccrImport}}.
#' @return Venn diagram plot. Can be plotted by \code{\link{grid.draw}}.
#' @export
#' @examples
#' data(hdacTR_smallExample)
#' trImported <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' vennPlot <- tppVenn(data=trImported)
#' grid.draw(vennPlot)

tppVenn <- function(data){
  if (is.data.frame(data)){
    data <- list(data)
  }
  
  ## 1.) Determine experiment names:
  expInfo <- sapply(data, annotation)
  expNames    <- expInfo["name",]
  names(data) <- expNames
  
  ## 2.) Determine experimental conditions and replicates (if available) 
  ## -> necessary to produce matched plot colors between Vehicle and Treatment
  grConditions <- expInfo["condition",]
  grReplicates <- as.numeric(expInfo["replicate",])
  
  ## 3.) Retrieve list of protein names per experiment:
  idsPerExperiment <- lapply(data, featureNames) # Holgers code: 'plotVennData'
  
  ## 4.) Determine plot colors, taking information about treatment group and 
  ## replicate into account (if available):
  plotCols <- plotColors(expConditions=grConditions,expReplicates=grReplicates)
  
  ## 5.) Fill annotation table:
  tableDF <- data.frame()
  for(en in expNames){
    lmnts <- featureNames(data[[en]])
    row_df = data.frame(condition=en, proteins=(length(lmnts)))
    tableDF = rbind(tableDF, row_df)
  }
  
  ## 7.) Determine plot header:
  mainTxt = ""
  
  ## II.) Produce plot:
  plotObj <- venn.diagram(x  = idsPerExperiment, 
                          filename = NULL, 
                          fill     = plotCols,
                          main            = mainTxt,
                          main.pos        = c(0.5, 0.05),
                          main.cex        = 1.5,
                          main.fontfamily = "sans",
                          lty   = 1,
                          lwd   = 0.25,
                          alpha = 0.3,
                          cex   = 2,
                          label.col        = 'blue',
                          label.fontfamily = "sans",
                          cat.cex        = 1.5,
                          cat.pos        = 0,
                          cat.dist       = 0.1,
                          cat.fontface   = 1,
                          cat.fontfamily = "sans",
                          cat.col        = plotCols) 
  
  ## III.) Add table that summarizes the proteins per group:
  plotObj = addTableToPlot(plotObj  = gTree(children=plotObj), 
                           tableDF  = tableDF, 
                           clrs     = plotCols)  
  return(plotObj)
}