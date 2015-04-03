#' Fit melting curves to all proteins in a dataset.
#' @details If the melting curve fitting procedure does not converge, it will be
#'   repeatedly started from perturbed starting parameters (maximum iterations 
#'   defined by argment \code{maxAttempts})
#'   
#' @return A list of ExpressionSets storing the data together with the melting
#'   curve parameters for each treatment condition and biological replicate.
#'   Each ExpressionSet contains the measured fold changes, as well as row and
#'   column metadata. In each ExpressionSet \code{S}, the fold changes can be
#'   accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. TMT labels and the corresponding temperatures are 
#'   returned by \code{S$labels} and \code{S$temperatures}.

#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=tpptrDefaultNormReqs())
#' normalizedData <- tpptrNorm$normData
#' hdacSubsets <- lapply(normalizedData, 
#'                       function(d) d[grepl("HDAC", featureNames(d))])
#' tpptrFittedHDACs <- tpptrCurveFit(hdacSubsets)
#' 
#' @param data list of \code{ExpressionSet}s with protein fold changes for curve fitting.
#' @param resultPath location where to store the melting curve plots.
#' @param ggplotTheme ggplot theme for melting curve plots.
#' @param doPlot boolan value indicating whether melting curves should 
#' be plotted, or whether just the curve parameters should be returned.
#' @param startPars start values for the melting curve parameters. Will be passed
#' to function \code{\link{nls}} for curve fitting.
#' @param maxAttempts maximal number of curve fitting attempts if model does not converge.
#' @param nCores either a numerical value given the desired number of CPUs, or 'max'
#' to automatically assign the maximum possible number (default).
#' @details The melting curve plots will be stored in a subfolder with name 
#' \code{Melting_Curves} at the location specified by \code{resultPath}. 
#' @seealso \code{\link{tppDefaultTheme}}
#' @export
tpptrCurveFit <- function(data, resultPath=NULL, ggplotTheme=tppDefaultTheme(), 
                          doPlot=TRUE, startPars=c("Pl"=0, "a"=550, "b"=10), 
                          maxAttempts=500, nCores='max'){
  if (is.null(resultPath)){
    doPlot <- FALSE
  }
  ## Check if output directory exists already. If not, create it here.
  if (doPlot){
    plotDir <- "Melting_Curves"
    if (!file.exists(file.path(resultPath, plotDir))) 
      dir.create(file.path(resultPath, plotDir), recursive=TRUE)    
  }
  
  ## Extract metadata about experiment names, conditions, replicates, proteins
  expInfo <- sapply(data, function(set){set@annotation})
  expNames        <- expInfo["name",]
  names(data) <- expNames
  grConditions <- expInfo["condition",]
  grReplicates <- as.numeric(expInfo["replicate",])
  protIDs <- unique(unlist(lapply(data, featureNames)))
  
  ## Prepare file names for melting curve plots. Replace special characters by
  ## '_'.
  if (doPlot){
    plotPathsFull <- file.path(plotDir, paste("meltCurve_",
                                              gsub("([^[:alnum:]])", "_", 
                                                   protIDs),".pdf", sep=""))
  } else {
    plotPathsFull <- rep("", length(protIDs))
  }
  names(plotPathsFull) <- protIDs
  
  ## Extract fold change and temperature values per protein and group:
  xMat <- t(sapply(data, function(set){set$temperature}))
  yMat <- data.frame(matrix(nrow=0, ncol=(ncol(data[[1]])+2)))
  for (g in expNames){
    yTmp <- exprs(data[[g]])
    if (nrow(yTmp)==0){
      yTmp <- matrix(NA_real_, nrow=length(protIDs), ncol=ncol(data[[g]]), 
                     dimnames=list(protIDs, colnames(yTmp)))
    }
    yMat <- rbind(yMat, data.frame("expName"=g, "protID"=rownames(yTmp), 
                                   "FC"=yTmp, stringsAsFactors=FALSE))
  }
  ## Determine operating system to determine whether legends should be added:
  osType <- Sys.info()['sysname']
  if (osType %in% c("Linux", "Windows")){
    addLegend <- TRUE    
  } else addLegend <- FALSE
  
  ## Calculate melting curves and plot results:
  message("Fitting melting curves to ", length(protIDs), " proteins.")
  nCores <- checkCPUs(cpus=nCores)
  doParallel::registerDoParallel(cores=nCores)
  t1 <- Sys.time()
  dfCurvePars <- foreach(p=protIDs, .combine=rbind, .inorder=FALSE, .verbose=FALSE, 
                         .packages = c("ggplot2", "gridExtra")) %dopar%
    fitMeltCurvesToProtein(xMat, yDF=subset(yMat, protID==p), startPars=startPars,
               maxAttempts=maxAttempts, expNames=expNames, resultPath=resultPath, 
               plotPathRel=plotPathsFull[[p]], protID=p, plotTheme=ggplotTheme, 
               grConds=grConditions, grReps=grReplicates, doPlot=doPlot, addLegend=addLegend)
  stopImplicitCluster()
  timeDiff <- Sys.time()-t1
  message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  gc(verbose=FALSE)
  message("Melting curves fitted sucessfully!")  
  
  ## Inform user about success rate:
  conv <- dfCurvePars$model_converged
  expr <- dfCurvePars$sufficient_data_for_fit
  message(sum(conv, na.rm=TRUE), " out of ", sum(expr), 
          " models with sufficient data points converged (",
          round(sum(conv, na.rm=TRUE)/sum(expr) * 100, 2)," %).\n")
  
  ## Store curve parameters in expression set annotation:
  data <- storeMeltCurveParams(data=data, params=dfCurvePars)  
  return(data)
}

stopImplicitCluster <- function()
{
    .options <- doParallel:::.options 
    if(exists(".revoDoParCluster", where=.options) && 
      !is.null(.options[['.revoDoParCluster']]))
    {
            stopCluster(.options[['.revoDoParCluster']])
            remove('.revoDoParCluster', envir=.options)
    }
}

