#' @title Plot dose response curves
#' @description \code{tppccrPlotCurves} plots the logistic dose response curves, 
#' as well as the underlying fold
#'   change measurements for each TPP-CCR experiment in a study.
#' @param data list of expressionSet objects containing protein fold changes, as 
#' well as fitted curve parameters.
#' @param fcTable optional long table with fold changes for each experiment. 
#' Can be provided instead of the input argument \code{data}.
#' @param curvePars optional long table of curve parameters per protein and 
#' experiment. Can be provided instead of the input argument \code{data}.
#' @param resultPath location where to store dose-response curve plots.
#' @param ggplotTheme ggplot theme for dose response curve plots.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param verbose print name of each plotted protein to the command line as a 
#' means of progress report.
#' 
#' @details \code{data} is a list of expressionSet objects created by
#' \code{\link{tppccrCurveFit}}. It contains 
#' the isobaric labels and administered drug concentrations in the 
#' \code{phenoData} and user-defined protein properties (including dose response 
#' curve parameters) in the \code{featureData}. Protein IDs are stored in the 
#' \code{featureNames}.
#' 
#' Measurements and compound effects for curve fitting can be provided 
#' by the arguments \code{fcTable} and \code{cpdEffects}, instead of being 
#' stored in expressionSets in \code{data}. 
#'   
#' If specified, \code{fcTable} needs to be a long 
#' table with column names "id" (the protein names), "concentration" (the fold 
#' changes), "labelName" (the isobaric label to each measurement), and 
#' "experiment" (e.g. "Vehicle_1" or "Panobinostat_1").
#'   
#' If specified, \code{curvePars} needs to be a long 
#' table with column names "id" (the protein names), "param" (curve parameter 
#' per protein and experiment, see TPP:::drCurveParamNames(names=TRUE, 
#' info=FALSE) for possibilities), and 
#' "experiment" (e.g. "Vehicle_1" or "Panobinostat_1").
#' 
#' The dose response curve plots will be stored in a subfolder with name 
#' \code{DoseResponse_Curves} at the location specified by \code{resultPath}.
#' 
#' @return A list of expressionSet objects storing fold changes,
#' as well as row and column metadata. In each expressionSet \code{S}, the fold 
#' changes
#'   can be accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. Isobaric labels and the corresponding 
#'   concentrations are 
#'   returned by \code{S$label} and \code{S$concentration}. Paths to the 
#'   produced plots are stored in code{featureData(S)$plot}.
#'  
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, 
#'                            data=hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' tppccrTransformed <- tppccrTransform(data=tppccrNorm)
#' tppccrFitted <- tppccrCurveFit(data=tppccrTransformed, nCores=1)
#' hdacSubset <- sapply(tppccrFitted, function(d)d[grepl("HDAC", rownames(d)),])
#' tppccrPlotted <- tppccrPlotCurves(hdacSubset, resultPath=getwd(), nCores = 1)
#' 
#' @seealso \code{\link{tppccrCurveFit}},\code{\link{tppDefaultTheme}}
#' 
#' @export
tppccrPlotCurves <- function(data=NULL, fcTable=NULL, curvePars=NULL,
                             resultPath=NULL, ggplotTheme=tppDefaultTheme(),
                             nCores="max",  verbose=FALSE){
  ## 1. Check if output directory exists already. If not, create it here.
  plotDir <- "DoseResponse_Curves"
  if (!file.exists(file.path(resultPath, plotDir))) {
    dir.create(file.path(resultPath, plotDir), recursive=TRUE)
  }
  
  ## 2. Define plot theme:
  theme_set(ggplotTheme)
  
  ## 3. Obtain long tables with fold changes and curve parameters
  if (!is.null(data)){
    isESetList <- ifelse (is.list(data)&identical(unique(sapply(data, class)),
                                                  "ExpressionSet"), TRUE, FALSE)
    if (isESetList) {
      fcTable <- eSetsToLongTable_fc(data)
      colnames(fcTable)[grep("labelValue",colnames(fcTable))] <- "concentration"
      
      fDatTable <- eSetsToLongTable_fData(data)
      parNames <- drCurveParamNames(names=TRUE, info=FALSE)
      curvePars <- subset(fDatTable, variable %in% parNames)
      curvePars$value <- as.numeric(curvePars$value)
      colnames(curvePars)[grep("variable", colnames(curvePars))] <- "param"
    }
  } else if (is.null(fcTable) | is.null(curvePars)) {
    stop("Please specify either 'data', or both 'fcTable' and 'curvePars'.")
  }
  
  ## 2. Ignore proteins with NA fold changes only
  idsValid <- fcTable %>% select(id, foldChange) %>% na.omit %>% 
    extract2("id") %>% unique %>% as.character
  fcFiltered <- subset(fcTable, id %in% idsValid) %>% 
    mutate(id = factor(as.character(id)))
  
  ## 3. Start parallelized DR plotting over all proteins:
  fcSplit  <- split(fcFiltered, fcFiltered$id)
  parSplit <- split(curvePars, curvePars$id)
  expNames <- unique(fcFiltered$experiment)
  message("Plotting dose response curves for ", length(fcSplit), " proteins.")  
  
  ## Determine operating system to decide whether legends should be added (on 
  ## Mac OS, legends interfere with parallelization due to a bug in one of the 
  ## involved packages):
  addLegend <- checkIfLegendPossible()
  
  ## Determine plot colors:
  plotCols <- plotColors(expConditions=c(NA,NA),  comparisonNums=c(NA,NA))
  
  
  nCores <- checkCPUs(cpus=nCores)
  t1 <- Sys.time()
  if (nCores == 1){
    plotFileNames <- foreach(pID=names(fcSplit), .combine=rbind, .inorder=FALSE, 
                             .verbose=FALSE)  %do% {
                               fcDF    = fcSplit[[pID]]
                               parDF   = parSplit[[pID]]
                               plotDRCurve(protID  = pID,
                                           fcDF    = fcDF,
                                           parDF   = parDF,
                                           plotDir = file.path(resultPath, 
                                                               plotDir),
                                           allExp  = expNames,
                                           addLegend = addLegend,
                                           plotCols = plotCols,
                                           verbose = verbose)
                             }
  } else if (nCores > 1){
    doParallel::registerDoParallel(cores=nCores)
    plotFileNames <- foreach(pID=names(fcSplit), .combine=rbind, .inorder=FALSE, 
                             .verbose=FALSE)  %dopar% {
                               fcDF    = fcSplit[[pID]]
                               parDF   = parSplit[[pID]]
                               plotDRCurve(protID  = pID,
                                           fcDF    = fcDF,
                                           parDF   = parDF,
                                           plotDir = file.path(resultPath, 
                                                               plotDir),
                                           allExp  = expNames,
                                           addLegend = addLegend, 
                                           plotCols = plotCols,
                                           verbose = verbose)
                             }
    stopImplicitCluster()
  }
  timeDiff <- Sys.time()-t1
  message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 2), " ", 
          units(timeDiff), "\n")
  gc(verbose=FALSE)
  message("Dose response curves plotted sucessfully!")  
  
  ## Store plot paths in the featureData of each expressionSet:
  plotFileNames$path <- file.path(plotDir, plotFileNames$path)
  if (!is.null(data)){
    for (e in expNames){
      datTmp <- data[[e]]
      plotDF <- join(data.frame("Protein_ID"=featureNames(datTmp)), 
                     plotFileNames, 
                     by="Protein_ID")
      featureData(data[[e]])$plot <- as.character(plotDF$path)
    }
    return(data)    
  } else {
    return(plotFileNames)
  }
}