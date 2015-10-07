#' @title Fit dose response curves
#' @description \code{tppccrCurveFit} fits logistic dose response curves to fold
#'   change measurements of a TPP-CCR experiment.
#' @param data list of expressionSet objects containing protein fold changes for 
#' dose response curve fitting.
#' @param fcTable optional long table with fold changes for each experiment. 
#' Can be provided instead of the input argument \code{data}.
#' @param cpdEffects optional long table of compound effects per protein and 
#' experiment. Can be provided instead of the input argument \code{data}.
#' @param slopeBounds bounds on the slope parameter for dose response curve
#'   fitting.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param verbose print name of each fitted protein to the command 
#' line as a means of progress report.
#' 
#' @details \code{data} is a list of expressionSet objects created by 
#'   \code{\link{tppccrImport}}. If desired, it can be already preprocessed by 
#'   \code{\link{tppccrNormalize}} or \code{\link{tppccrTransform}}. It contains
#'   the isobaric labels and administered drug concentrations in the 
#'   \code{phenoData} and user-defined protein properties in the 
#'   \code{featureData}. Protein IDs are stored in the \code{featureNames}.
#'   
#'   Measurements and compound effects for curve fitting can be provided 
#'   by the arguments \code{fcTable} and \code{cpdEffects}, instead of being 
#'   stored in expressionSets in \code{data}. 
#'   
#'   If specified, \code{fcTable} needs to be a long 
#'   table with column names "id" (the protein names), "concentration" (the fold 
#'   changes), "labelName" (the isobaric label to each measurement), and 
#'   "experiment" (e.g. "Vehicle_1" or "Panobinostat_1").
#'   
#'   If specified, \code{cpdEffects} needs to be a long 
#'   table with column names "id" (the protein names), "cpdEff" (character 
#'   vector of compound effects, may contain NAs), and 
#'   "experiment" (e.g. "Vehicle_1" or "Panobinostat_1").
#'   
#' @return A list of expressionSet objects storing fold changes, the fitted
#'   curve parameters, as well as row and column metadata. In each expressionSet
#'   \code{S}, the fold changes can be accessed by \code{exprs(S)}. Protein
#'   expNames can be accessed by \code{featureNames(S)}. Isobaric labels and the
#'   corresponding concentrations are returned by \code{S$label} and
#'   \code{S$concentration}. The fitted curve parameters are stored in
#'   code{featureData(S)}.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, 
#'                            data=hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' tppccrTransformed <- tppccrTransform(data=tppccrNorm)
#' tppccrFitted <- tppccrCurveFit(data=tppccrTransformed, nCores=1)
#' 
#' @seealso \code{\link{tppccrImport}}, \code{\link{tppccrNormalize}}, 
#'   \code{\link{tppccrTransform}}
#'   
#' @export
tppccrCurveFit <- function(data=NULL, fcTable=NULL, cpdEffects=NULL, 
                           slopeBounds=c(1,50), nCores='max', verbose=FALSE){
  
  ## 1. Obtain long tables with fold changes and compound effects
  if (!is.null(data)){
    isESetList <- ifelse (is.list(data)&identical(unique(sapply(data, class)),
                                                  "ExpressionSet"), TRUE, FALSE)
    if (isESetList) {
      fcTable <- eSetsToLongTable_fc(data)
      colnames(fcTable)[grep("labelValue", colnames(fcTable))] <- "concentration"
      
      fDatTable <- eSetsToLongTable_fData(data)
      cpdEffects <- subset(fDatTable, variable=="compound_effect", 
                           select=c("id", "value", "experiment"))
      colnames(cpdEffects)[grep("value", colnames(cpdEffects))] <- "cpdEff"
    }
  } else if (is.null(fcTable) | is.null(cpdEffects)) {
    stop("Please specify either 'data', or both 'fcTable' and 'cpdEffects'.")
  }
  fcTable$experiment <- as.character(fcTable$experiment)
  fcTable$id         <- as.character(fcTable$id)
  fcTable$uniqueID    <- paste(fcTable$id, fcTable$experiment, sep="_")
  cpdEffects$uniqueID <- paste(cpdEffects$id, cpdEffects$experiment, sep="_")
  
  ## 2. Ignore proteins with NA fold changes only
  fcNonNA <- c()
  for (en in unique(fcTable$experiment)){
    fcTmp <- subset(fcTable, experiment==en)
    idsValid <- unique(fcTmp[which(!is.na(fcTmp$foldChange)),]$id)
    fcFiltered <- subset(fcTmp, id %in% idsValid)
    fcNonNA <- rbind(fcNonNA, fcFiltered)
  }
  
  ## 3. Start parallelized DR curve fitting over all proteins and experiments:
  fcSplit  <- split(fcNonNA, fcNonNA$uniqueID)
  effSplit <- split(cpdEffects, cpdEffects$uniqueID)
  
  message("Fitting ",length(fcSplit), " individual dose response curves to ", 
          length(unique(fcNonNA$id))," proteins.")  
  nCores <- checkCPUs(cpus=nCores)
  t1 <- Sys.time()
  if (nCores == 1){
    parsFitted <- foreach(i=names(fcSplit), .combine=rbind, .inorder=FALSE, 
                          .verbose=FALSE)  %do%
      fitDRCurve(protID   = unique(fcSplit[[i]]$id), 
                 expName  = unique(fcSplit[[i]]$experiment), 
                 dose     = fcSplit[[i]]$concentration, 
                 response = fcSplit[[i]]$foldChange, 
                 cpd_effect = effSplit[[i]]$cpdEff,
                 slBds     = slopeBounds, 
                 verbose   = verbose)
  } else if (nCores > 1){
    doParallel::registerDoParallel(cores=nCores)
    parsFitted <- foreach(i=names(fcSplit), .combine=rbind, .inorder=FALSE, 
                          .verbose=FALSE)  %dopar%
      fitDRCurve(protID   = unique(fcSplit[[i]]$id), 
                 expName  = unique(fcSplit[[i]]$experiment), 
                 dose     = fcSplit[[i]]$concentration, 
                 response = fcSplit[[i]]$foldChange, 
                 cpd_effect = effSplit[[i]]$cpdEff,
                 slBds   = slopeBounds, 
                 verbose  = verbose)
    stopImplicitCluster()    
  }
  timeDiff <- Sys.time()-t1
  message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 2), " ", 
          units(timeDiff), "\n")
  gc(verbose=FALSE)
  message("Dose response curves fitted sucessfully!")  
  
  ## 4. Inform user about success rate:
  conv <- parsFitted$model_converged
  expr <- parsFitted$sufficient_data_for_fit
  message(sum(conv, na.rm=TRUE), " out of ", sum(expr), 
          " models with sufficient data points converged (",
          round(sum(conv, na.rm=TRUE)/sum(expr) * 100, 2)," %).\n")
  
  ## 5. Store curve parameters in featureData (will be used to re-compute the 
  ##    curves for plotting)
  if (!is.null(data)){
    dataFitted <- storeDRCurveParams(data=data, params=parsFitted)
    return(dataFitted)    
  } else {
    return(parsFitted)
  }
}