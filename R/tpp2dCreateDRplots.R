#' @title Create dose response curve plots for 2D-TPP data
#' @description Generates a list of dose response curve plots per protein and 
#' temperature point.
#'  
#' @return A list of successfully generated plot objects of class 
#'         \code{'ggplot'}
#' @examples
#' data(panobinostat_2DTPP_smallExample)
#' 
#' # Import data:
#' datIn <- tpp2dImport(configTable = panobinostat_2DTPP_config,
#'                       data = panobinostat_2DTPP_data,
#'                       idVar = "representative",
#'                       addCol = "clustername",
#'                       intensityStr = "sumionarea_protein_",
#'                       nonZeroCols = "qusm")
#' 
#' # Compute fold changes:
#' fcData2d <- tpp2dComputeFoldChanges(data = datIn)
#' normData2d <- tpp2dNormalize(data = fcData2d)
#' ccr2dResults <- tpp2dCurveFit(data = normData2d)
#' allCurves <- tpp2dCreateDRplots(data = ccr2dResults, type = "all")
#' allCurves[["HDAC1"]]
#' 
#' 
#' @param data the data that should be plotted.
#' @param verbose boolean variable stating whether a print description of 
#' problems/success for plotting of each protein should be printed.
#' @param type string defining which curves to display (see details).
#' @param paletteName color palette (see details).
#' 
#' @details 
#' \code{data} is a data frame in wide table format returned by function 
#' \code{\link{tpp2dCurveFit}}. Its attributes contain information about the 
#' experiment names, temperatures, isobaric labels, as well as instructions on 
#' how to find the relevant columns in the wide table.
#' 
#' \code{type} defines which curves to display per plot. Possible values are:
#'   \itemize{ 
#'   \item{"all": Create one plot per protein. This plot simultaneously 
#'   displays the curves for all available temperatures for this protein 
#'   (the default).} 
#'   \item{"good": Create one plot per protein. This plot displays all
#'   dose response curves with a high goodness-of-fit.
#'   Choose this option to save runtime by focusing only on the reliable fits.} 
#'   \item{"single": Create one separate plot per protein and temperature. 
#'   This plot displays all dose response curves with a high goodness-of-fit.} 
#'   }
#'  
#'  \code{paletteName} specifies the color palette to be used by the \code{\link{brewer.pal}} 
#' function from the \code{RColorBrewer} package to assign a separate color to 
#' each concentration.
#'   
#' @seealso \code{\link{tpp2dCurveFit}} \code{\link{brewer.pal}}
#' @export

tpp2dCreateDRplots <- function(data = NULL, type = "all", verbose = FALSE, 
                               paletteName = "Spectral"){
  # problem: update documentation of 'good' option.
  # suggestions: 
  # for which the \eqn{R^2} value of the fit exceeded 
  #   \eqn{0.8}.
  #   (Corresponds to the \code{passed_filter} column in the output of 
  #   \code{\link{tpp2dCurveFit}}).
  
  # Obtain config table used for data import (stored as attribute of imported data):
  configTable <- attr(data, "configTable")
  
  # Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(data, "importSettings")
  proteinIdCol <- importSettings$proteinIdCol
  
  # Choose correct fold change column prefix (automatically detects whether
  # to use the prefix for normalized columns).
  finalFcPrefix <- obtain_fcStr_from_df_annotation(dat = data)
  
  if (type == "all"){
    # generate joint plots for all proteins detected
    plotList <- helperFctPlotAll(configTable = configTable, 
                                 dataTable = data, 
                                 idVar = proteinIdCol, 
                                 fcStr = finalFcPrefix,
                                 verbose = verbose,
                                 paletteName = paletteName)
    
  } else if (type == "single"){
    # generate single plots for all protein in each condition fitted with sufficient R2
    plotList <- helperFctPlotSingle(configTable = configTable, 
                                    dataTable = data, 
                                    idVar = proteinIdCol, 
                                    fcStr = finalFcPrefix,
                                    verbose = verbose)
    
  } else if (type == "good"){
    # generate joint plots for all proteins detected with sufficient R2
    plotList <- helperFctPlotGood(configTable = configTable, 
                                  dataTable = data, 
                                  idVar = proteinIdCol, 
                                  fcStr = finalFcPrefix,
                                  verbose = verbose,
                                  paletteName = paletteName)
    
  } else {
    warning("No plots could be produced because of invalid 'type' argument.",
            "Given value: 'type' = '", type, "'.")
    plotList <- list()
  }
  return(plotList)
}
