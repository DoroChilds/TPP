#' Fit melting curves to all proteins in a dataset.
#' @details If the melting curve fitting procedure does not converge, it will be
#'   repeatedly started from perturbed starting parameters (maximum iterations 
#'   defined by argument \code{maxAttempts})
#'   
#'   If \code{doPlot = TRUE}, melting curves are be plotted in individual files
#'   per protein. Each file is named by its unique identifier. Filenames are
#'   truncated to 255 characters (requirement by most operation systems). 
#'   Truncated filenames are indicated by the suffix "_truncated[d]", where [d] 
#'   is a unique number to avoid redundancies.
#'   
#' @return A list of ExpressionSets storing the data together with the melting
#'   curve parameters for each experiment.
#'   Each ExpressionSet contains the measured fold changes, as well as row and
#'   column metadata. In each ExpressionSet \code{S}, the fold changes can be
#'   accessed by \code{Biobase::exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. Isobaric labels and the corresponding temperatures are 
#'   returned by \code{S$label} and \code{S$temperature}.

#' 
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=tpptrDefaultNormReqs())
#' normalizedData <- tpptrNorm$normData
#' hdacSubsets <- lapply(normalizedData, 
#'                       function(d) d[grepl("HDAC", Biobase::featureNames(d))])
#' tpptrFittedHDACs <- tpptrCurveFit(hdacSubsets, nCores=1)
#' # Show estimated parameters for vehicle and treatment experiments:
#' Biobase::pData(Biobase::featureData(tpptrFittedHDACs[["Vehicle_1"]]))
#' Biobase::pData(Biobase::featureData(tpptrFittedHDACs[["Panobinostat_1"]]))
#' 
#' @param data list of \code{ExpressionSet}s with protein fold changes for curve
#'   fitting.
#' @param dataCI list of \code{ExpressionSet}s with protein fold change confidence 
#' intervals for curve fitting. Default to NULL.
#' @param resultPath location where to store the melting curve plots.
#' @param ggplotTheme ggplot theme for melting curve plots.
#' @param doPlot boolean value indicating whether melting curves should be
#'   plotted, or whether just the curve parameters should be returned.
#' @param startPars start values for the melting curve parameters. Will be
#'   passed to function \code{\link{nls}} for curve fitting.
#' @param maxAttempts maximal number of curve fitting attempts if model does not
#'   converge.
#' @param nCores either a numerical value given the desired number of CPUs, or
#'   'max' to automatically assign the maximum possible number (default).
#' @param verbose plot name of each fitted protein to the command lin as a means
#'   of progress report.
#' @details The melting curve plots will be stored in a subfolder with name 
#'   \code{Melting_Curves} at the location specified by \code{resultPath}.
#' @seealso \code{\link{tppDefaultTheme}}
#' @export
tpptrCurveFit <- function(data, dataCI=NULL, resultPath=NULL, 
                          ggplotTheme=tppDefaultTheme(), doPlot=TRUE, 
                          startPars=c("Pl"=0, "a"=550, "b"=10), 
                          maxAttempts=500, nCores='max', verbose=FALSE){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  p <- NULL
  
  # set parameter if curves are plotted
  #doPlot <- getOption("TPPTR_plot") && !is.null(resultPath)
  doPlot <- doPlot && !is.null(resultPath)
  
  # get CI option: should CIs be used in the fit or not ?
  useCI <- getOption("TPPTR_CI")
  if(is.null(useCI)) useCI <- FALSE

  
  ## Extract metadata about experiment names, conditions, and proteins
  expInfo      <- sapply(data, annotation)
  expNames     <- expInfo["name",]
  grConditions <- expInfo["condition",]
  compDF       <- createComparisonTable(infoTable=expInfo)
  protIDs <- unique(unlist(lapply(data, featureNames)))
  
  ## Prepare plots (if doPlot=TRUE)
  if (doPlot){
    ## Check if output directory exists already. If not, create it here.
    plotDir <- "Melting_Curves"
    if (!file.exists(file.path(resultPath, plotDir))){
      dir.create(file.path(resultPath, plotDir), recursive=TRUE)
    }
    ## File names for melting curve plots: Replace special characters by '_'.
    fNames <- paste0("meltCurve_",gsub("([^[:alnum:]])","_", protIDs))

    # limit file name length to 255 characters to avoid crashes on most filesystems:
    maxLen <- 255 - nchar(".pdf") - nchar("_truncated") - nchar(as.character(length(fNames)))
    tooLong <- nchar(fNames) > maxLen
    cropSuffix <- paste0("_truncated", 1:sum(tooLong))
    fNames <- sapply(fNames, function(fTmp) {
      fNew <- substr(fTmp, 1, min(maxLen, nchar(fTmp)))
    }, simplify = TRUE, USE.NAMES = FALSE)
    fNames[tooLong] <- paste0(fNames[tooLong], cropSuffix)
    fNames <- paste0(fNames, ".pdf")
    
    
    plotPathsFull <- file.path(plotDir, fNames)
  } else {
    plotPathsFull <- rep("", length(protIDs))
  }
  names(plotPathsFull) <- protIDs
  
  ## Extract fold change and temperature values per protein and group:
  xMat <- t(sapply(data, function(set){set$temperature}))
  yMat <- data.frame(matrix(nrow=0, ncol=(ncol(data[[1]])+2)))
  for (g in expNames){
    yTmp <- Biobase::exprs(data[[g]])
    if (nrow(yTmp)==0){
      yTmp <- matrix(NA_real_, nrow=length(protIDs), ncol=ncol(data[[g]]), 
                     dimnames=list(protIDs, colnames(yTmp)))
    }
    ## Ensure that fold changes stay sorted by temperature and are not re-sorted
    ## by labels during the rbind command (important if different experiments 
    ## assign different temperature order to tmt labels)
    colnames(yTmp) <- 1:ncol(yTmp) 
    yMat <- rbind(yMat, data.frame("expName"=g, "protID"=rownames(yTmp), 
                                   "FC"=yTmp, stringsAsFactors=FALSE))
    
  }
  
  if (useCI){
    ciMat <- data.frame(matrix(nrow=0, ncol=(ncol(data[[1]])+2)))
    for (g in expNames){
      ciTmp <- Biobase::exprs(dataCI[[g]])
      if (nrow(ciTmp)==0){
        ciTmp <- matrix(NA_real_, nrow=length(protIDs), ncol=ncol(dataCI[[g]]), 
                        dimnames=list(protIDs, colnames(ciTmp)))
      }
      ## Ensure that fold changes stay sorted by temperature and are not re-sorted
      ## by labels during the rbind command (important if different experiments 
      ## assign different temperature order to tmt labels)
      colnames(ciTmp) <- 1:ncol(ciTmp) 
      ciMat <- rbind(ciMat, data.frame("expName"=g, "protID"=rownames(ciTmp), 
                                       "CI"=ciTmp, stringsAsFactors=FALSE))
      
    }
  }
  
  
  ## Determine operating system to decide whether legends should be added 
  ## (legends interfere with parallelization due to a bug in one of the involved packages):
  addLegend <- checkIfLegendPossible()
  
  ## Calculate melting curves and plot results:
  message("Fitting melting curves to ", length(protIDs), " proteins.")
  nCores <- checkCPUs(cpus=nCores)
  t1 <- Sys.time()
  
  
  # ciDF = ifelse(is.null(ciMat), NULL,  subset(ciMat, protID==p))
  
  # differentiate between parallel and serial execution:
  # if nCores == 1 use serial execution via lapply
  # if nCores > 1 use the foreach package for multicore execution
  if (nCores > 1) {
    doParallel::registerDoParallel(cores = nCores)
    
    # fitMeltCurves <- function(xMat, yDF, startPars, maxAttempts, expNames, 
    # protID, verbose)
    
    
    # parallel execution with foreach and dopar 
    dfCurvePars <- foreach(p=protIDs, .combine=rbind, .inorder=FALSE, .verbose=FALSE, 
                           .packages = c("ggplot2", "gridExtra")) %dopar% {
                             tpptrHelperFitandPlot(p, yMat, xMat, ciMat, 
                                                   startPars, maxAttempts, 
                                                   expNames, verbose, 
                                                   ggplotTheme, grConditions, 
                                                   compDF, 
                                                   addLegend, resultPath, 
                                                   plotPathsFull, useCI, doPlot)
                           }
    
    
    stopImplicitCluster()  
    # } <- deleted
    
  } else {
    
    # serial execution with laapply and do.call("rbind", ...) 
    dfCurvePars <- lapply(protIDs, 
                          function(p){          
                            tpptrHelperFitandPlot(p, yMat, xMat, ciMat, 
                                                  startPars, maxAttempts, 
                                                  expNames, verbose, 
                                                  ggplotTheme, grConditions, 
                                                  compDF, 
                                                  addLegend, resultPath, 
                                                  plotPathsFull, useCI, doPlot)
                          })
    
    dfCurvePars = do.call('rbind', dfCurvePars)
  }
  
  
  timeDiff <- Sys.time()-t1
  message("Runtime (", nCores, " CPUs used): ", round(timeDiff, 2), " ", 
          units(timeDiff), "\n")
  gc(verbose=FALSE)
  message("Melting curves fitted sucessfully!")  
  
  ## Inform user about success rate:
  dfCurvePars$sufficient_data_for_fit <- as.logical(dfCurvePars$sufficient_data_for_fit)
  dfCurvePars$model_converged         <- as.logical(dfCurvePars$model_converged)
  
  conv <- dfCurvePars$model_converged
  expr <- dfCurvePars$sufficient_data_for_fit
  message(sum(conv, na.rm=TRUE), " out of ", sum(expr), 
          " models with sufficient data points converged (",
          round(sum(conv, na.rm=TRUE)/sum(expr) * 100, 2)," %).\n")
  
  ## Store curve parameters in expression set annotation:
  data <- storeMeltCurveParams(data=data, params=dfCurvePars)  
  return(data)
}


tpptrHelperFitandPlot <- function(p, yMat, xMat, ciMat, startPars, maxAttempts, 
                                  expNames, verbose, ggplotTheme, grConditions, 
                                  compDF, addLegend, resultPath, plotPathsFull, 
                                  useCI, doPlot){
  # Helper function that combines the code for fitting and plotting of melting 
  # curves for parallel or sequential (in lapply) execution
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  protID <- NULL
  
  # retrieving options via getOptions doesn't work with parallel execution
  # therefore options are passed as variables (useCI, doPlot)
  
  yDF = subset(yMat, protID==p)
  
  resFC <- fitMeltCurves(xMat, yDF=yDF, colPrefix = "FC",
                         startPars=startPars,maxAttempts=maxAttempts, 
                         expNames=expNames, 
                         protID=p, verbose=verbose)
  
  # if we ant to use the CIs give the subseit of the CI matrix to the function
  # otherwise NULL
  if (useCI) {
    ciDF=subset(ciMat, protID==p)
    ciDF=ciDF[match(yDF$expName, ciDF$expName),]
    
    ciDFUpper <- ciDFLower <- ciDF
    ciDFUpper[,-(1:2)] <- yDF[,-(1:2)] + ciDF[,-(1:2)] / 2
    ciDFLower[,-(1:2)] <- yDF[,-(1:2)] - ciDF[,-(1:2)] / 2
    
    resUpper <- fitMeltCurves(xMat, yDF=ciDFUpper, colPrefix = "CI",
                              startPars=startPars,maxAttempts=maxAttempts, 
                              expNames=expNames, 
                              protID=p, verbose=verbose)
    
    resLower <- fitMeltCurves(xMat, yDF=ciDFLower, colPrefix = "CI",
                              startPars=startPars,maxAttempts=maxAttempts, 
                              expNames=expNames, 
                              protID=p, verbose=verbose)
    
    listUpper = resUpper[[3]]
    listLower = resLower[[3]]
    
    CI_reportData = data.frame(CI_meltPointUpper = resUpper[[1]]$meltPoint,
                               CI_meltPointLower = resLower[[1]]$meltPoint,
                               CI_meltPoint_delta = resUpper[[1]]$meltPoint - resLower[[1]]$meltPoint)
    
  } else {  
    listUpper = listLower = NULL
  }
  
  if(doPlot){
    pl <- plotMeltingCurve(modelList = resFC[[3]], 
                           listUpper = listUpper, 
                           listLower = listLower,
                           xMat = xMat, 
                           fcMat = resFC[[2]], 
                           curvePars = resFC[[1]], 
                           protID = p,
                           plotTheme = ggplotTheme, 
                           expConditions = grConditions, 
                           expComps = compDF, 
                           addLegend = addLegend, 
                           useCI=useCI)
    
    if(is.null(pl)){
      plotPathRel <- NA_character_
    } else {
      plotPathRel <- plotPathsFull[p]
      pdf(file=file.path(resultPath, plotPathRel), width=7.87, height=9.84, 
          useDingbats=FALSE)
      grid.draw(pl)
      dev.off()
      
    }
  } else {
    plotPathRel <- NA_character_
  }
  
  curveParsWholeProt <- resFC[[1]]
  curveParsWholeProt$protID <- p
  curveParsWholeProt$plot   <- plotPathRel
  curveParsWholeProt$expName <- expNames  
  
  if (useCI){
    curveParsWholeProt = cbind(curveParsWholeProt, CI_reportData)
  }
  
  curveParsWholeProt
}

