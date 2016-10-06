#' @title Analyze TPP-TR experiment
#'   
#' @description Performs analysis of a TPP-TR experiment by invoking routines 
#'   for data import, data processing, normalization, curve fitting, and 
#'   production of the result table.
#'   
#' @references Savitski, M. M., Reinhard, F. B., Franken, H., Werner, T., 
#'   Savitski, M. F., Eberhard, D., ... & Drewes, G. (2014). Tracking cancer 
#'   drugs in living cells by thermal profiling of the proteome. Science, 
#'   346(6205), 1255784.
#'   
#' @return A data frame in which the fit results are stored row-wise for each 
#'   protein.
#'   
#' @details Invokes the following steps: \enumerate{ \item Import data using the
#'   \code{\link{tpptrImport}} function. \item Perform normalization (optional)
#'   using the \code{\link{tpptrNormalize}} function. To perform normalization, 
#'   set argument \code{normalize=TRUE}. The normalization will be filtered 
#'   according to the criteria specified in the \code{normReqs} argument (also 
#'   see the documentation of \code{\link{tpptrNormalize}} and 
#'   \code{\link{tpptrDefaultNormReqs}} for further information). \item Fit
#'   melting curves using the function \code{\link{tpptrCurveFit}}. \item
#'   Produce result table using the function \code{\link{tpptrAnalyzeMeltingCurves}}. 
#'   \item Export results to Excel using the function \code{\link{tppExport}}. }
#'   
#'   The default settings are tailored towards the output of the python package 
#'   isobarQuant, but can be customised to your own dataset by the arguments 
#'   \code{idVar, fcStr, naStrs, qualColName}.
#'   
#'   If \code{resultPath} is not specified, the location of the first input file
#'   specified in \code{configTable} will be used. If the input data are not 
#'   specified in \code{configTable}, no result path will be set. This means 
#'   that no output files or melting curve plots are produced and 
#'   \code{analyzeTPPTR} just returns the results as a data frame.
#'   
#'   The function \code{analyzeTPPTR} reports intermediate results to the 
#'   command line. To suppress this, use \code{\link{suppressMessages}}.
#'   
#'   The \code{configTable} argument is a dataframe, or the path to a 
#'   spreadsheet (tab-delimited text-file or xlsx format). Information about 
#'   each experiment is stored row-wise. It contains the following columns: 
#'   \itemize{ \item{\code{Path}:}{location of each datafile. Alternatively, 
#'   data can be directly handed over by the \code{data} argument.} 
#'   \item{\code{Experiment}: }{unique experiment names.} 
#'   \item{\code{Condition}: }{experimental conditions of each dataset.} 
#'   \item{Label columns: } each isobaric label names a column that contains the
#'   temperatures administered for the label in the individual experiments. }
#'   
#'   The argument \code{methods} can be one of the following:
#'   More than one method can be specified. For example, parametric testing of 
#'   melting points and nonparametric spline-based goodness-of-fit tests can be 
#'   performed seqeuentially in the same analysis. The results are then written 
#'   to separate columns of the output table.
#'   
#'   If \code{methods} contains "meltcurvefit", melting curve plots will be 
#'   stored in a subfolder with name \code{Melting_Curves} at the location 
#'   specified by \code{resultPath}.
#'   If \code{methods} contains "splinefit", plots of the natural spline fits will be 
#'   stored in a subfolder with name \code{Spline_Fits} at the location 
#'   specified by \code{resultPath}.
#'   
#'   The argument \code{nCores} could be either 'max' (use all available cores) 
#'   or an upper limit of CPUs to be used.
#'   
#'   If \code{doPlot = TRUE}, melting curve plots are generated seperately for 
#'   each protein and stored in separate pdfs.
#'   Each file is named by the unique protein identifier. Filenames are
#'   truncated to 255 characters (requirement by most operation systems). 
#'   Truncated filenames are indicated by the suffix "_truncated[d]", where [d] 
#'   is a unique number to avoid redundancies.
#'   All melting curve plots are stored in a subfolder with name 
#'   \code{Melting_Curves} at the location specified by \code{resultPath}.
#'   
#'   If the melting curve fitting procedure does not converge, it will be 
#'   repeatedly started from perturbed starting parameters (maximum iterations 
#'   defined by argument \code{maxAttempts}).
#'   
#'   Argument \code{splineDF} specifies the degrees of freedom for natural
#'   spline fitting. As a single numeric value, it is directly passed on to the
#'   \code{splineDF} argument of \code{splines::ns}. Experience shows that
#'   \code{splineDF = 4} yields good results for TPP data sets with 10
#'   temperature points. It is also possible to provide a numeric vector. In
#'   this case, splines are fitted for each entry and the optimal value is
#'   chosen per protein using Akaike's Information criterion.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrResults <- analyzeTPPTR(configTable = hdacTR_config, data = hdacTR_data, 
#'                 methods = "splinefit", nCores = 1)
#' 
#' @param configTable dataframe, or character object with the path to a file, 
#'   that specifies important details of the TPP-TR experiment. See Section 
#'   \code{details} for instructions how to create this object.
#' @param data single dataframe, or list of dataframes, containing fold change 
#'   measurements and additional annotation columns to be imported. Can be used 
#'   instead of specifying the file path in the \code{configTable} argument.
#' @param resultPath location where to store melting curve plots, intermediate 
#'   results, and the final results table.
#' @param methods statistical methods for modeling melting behaviour and detecting 
#'   significant differences between experimental conditions. Ich more than one 
#'   method are specified, results will be computed for each and concatenated in 
#'   the result talble (default: meltcurvefit).
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param ciStr character string indicating which columns contain confidence 
#' intervals for the fold change measurements. If specified, confidence 
#' intervals will be plotted around the melting curves.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param normalize perform normalization (default: TRUE).
#' @param normReqs list of filtering criteria for construction of the 
#'   normalization set.
#' @param ggplotTheme ggplot theme for melting curve plots.
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param startPars start values for the melting curve parameters. Will be 
#'   passed to function \code{\link{nls}} for curve fitting.
#' @param maxAttempts maximal number of curve fitting attempts if model does not
#'   converge.
#' @param plotCurves boolan value indicating whether melting curves should be 
#'   plotted. Deactivating plotting decreases runtime.
#' @param fixedReference name of a fixed reference experiment for normaliztion. 
#'   If NULL (default), the experiment with the best R2 when fitting a melting 
#'   curve through the median fold changes is chosen as the reference.
#' @param pValMethod Method for p-value computation. Currently restricted to 
#'   'robustZ' (see Cox & Mann (2008)).
#' @param pValParams optional list of parameters for p-value computation.
#' @param pValFilter optional list of filtering criteria to be applied before 
#'   p-value computation.
#' @param verbose print name of each fitted protein to the command lin as a
#'   means of progress report.
#' @param xlsxExport boolean value indicating whether to produce result table in
#'   .xlsx format (requires package \code{openxlsx} and a zip application to be 
#'   installed).
#' @param splineDF degrees of freedom for natural spline fitting.
#'   
#' @seealso tppDefaultTheme, tpptrImport, tpptrNormalize, tpptrCurveFit, 
#' tpptrAnalyzeMeltingCurves
#' @export
analyzeTPPTR <- function(configTable, data = NULL, resultPath = NULL, 
                         methods = c("meltcurvefit", "splinefit"),
                         idVar = "gene_name", fcStr = "rel_fc_", ciStr = NULL, # settings for data import
                         naStrs = c("NA", "n/d", "NaN", "<NA>"), 
                         qualColName = "qupm", 
                         normalize = TRUE, normReqs = tpptrDefaultNormReqs(), # settings for cross-experiment normalization
                         ggplotTheme = tppDefaultTheme(), 
                         nCores = 'max', 
                         startPars = c("Pl" = 0, "a" = 550, "b" = 10), # settings for melting curve fitting
                         splineDF = 4, # settings for spline fitting
                         maxAttempts = 500, plotCurves = TRUE, 
                         fixedReference = NULL, 
                         pValMethod = "robustZ", # to do: remove this argument
                         pValFilter = list(minR2 = 0.8, maxPlateau = 0.3), # settings for robust-z-test
                         pValParams = list(binWidth = 300), 
                         verbose = FALSE, xlsxExport = TRUE){
  message("This is TPP version ", packageVersion("TPP"),".")
  
  ## Import data:
  trData <- tpptrImport(configTable=configTable, data=data, idVar=idVar, fcStr=fcStr, 
                        naStrs=naStrs, qualColName=qualColName)
  
  if(!is.null(ciStr)){
    
    # set option to indiacte use of confidence intervals during analysis
    options("TPPTR_CI" = TRUE)
    
    trDataCI <- tpptrImport(configTable=configTable, data=data, idVar=idVar, fcStr=ciStr, 
                            naStrs=naStrs, qualColName=qualColName)
    
    stopifnot(all(sapply(names(trData), function(i){all(dim(trData[[i]]) == dim(trDataCI[[i]]))})))
    
  } else {
    options("TPPTR_CI" = FALSE)
    trDataCI <- NULL
  }
  
  # # set plot parameter for curve plotting
  # options("TPPTR_plot" = plotCurves)
  
  expInfo   <- sapply(trData, annotation)
  expNames  <- names(trData)
  expNum   <- length(expNames)
  
  expConds   <- sapply(trData, function(s) s@annotation[["condition"]])
  expComps  <- createComparisonTable(infoTable=expInfo)  
  
  ## Extract directory from the filenames in config table, if specified:
  confgFields <- suppressMessages(importCheckConfigTable(infoTable=configTable, 
                                                         type="TR"))
  files      <- confgFields$files
  outDirList <- importFct_makeOutputDirs(outDir=resultPath, fNames=files)
  flagDoWrite <- outDirList$doWrite
  pathDataObj <- outDirList$pathDataObj
  pathExcel=resultPath <- outDirList$outDir
  if (!flagDoWrite) plotCurves <- FALSE
  
  ## Save imported data before normalization:
  if (flagDoWrite){
    save(list=c("trData", "trDataCI"), file=file.path(pathDataObj, "importedData.RData"))
  }
  
  ## Normalize data:
  if (normalize){
    normResults <- tpptrNormalize(data=trData, normReqs=normReqs, 
                                  qcPlotTheme=ggplotTheme, qcPlotPath=NULL, fixedReference=fixedReference)
    trDataNormalized <- normResults[["normData"]]
    ## Save normalized data:
    if (flagDoWrite){
      save(list=c("trDataNormalized"), file=file.path(pathDataObj, 
                                                      "normalizedData.RData"))    
    }
  } else {
    trDataNormalized <- trData
  }
  
  allIDs <- lapply(trDataNormalized, featureNames) %>% unlist %>% unname %>% 
    unique %>% sort
  resultTable <- data.frame(Protein_ID = allIDs, stringsAsFactors = FALSE)
  
  if (any(methods == "meltcurvefit")){
    ## Fit melting curves:
    trDataFitted <- tpptrCurveFit(data=trDataNormalized, dataCI=trDataCI,
                                  resultPath=resultPath,
                                  ggplotTheme=ggplotTheme, doPlot=plotCurves,
                                  startPars=startPars, maxAttempts=maxAttempts, 
                                  nCores=nCores, verbose=verbose)
    ## Save data including parameters of the model fits:
    if (flagDoWrite){
      save(list=c("trDataFitted"), file=file.path(pathDataObj, "fittedData.RData"))
    }
    
    ## Analyse melting curves and create result table:
    meltCurveResultTable <- tpptrAnalyzeMeltingCurves(data = trDataFitted, 
                                                      pValMethod = pValMethod, 
                                                      pValFilter = pValFilter, 
                                                      pValParams = pValParams)
    if (flagDoWrite){
      save(list = c("meltCurveResultTable"), 
           file = file.path(pathDataObj, "trResultsMeltCurveFit.RData")) # "results_TPP_TR.RData"))
    }
    resultTable <- left_join(resultTable, meltCurveResultTable, by = "Protein_ID")
  }
  
  if (any(methods == "splinefit")){
    ## Preparation: convert eSets to long tables:
    trDataLongTables <- trDataNormalized %>% tpptrTidyUpESets
    annotData <- trDataLongTables %>% extract2("proteinAnnotation")
    fitData <- trDataLongTables %>% extract2("proteinMeasurements")
    
    splineFitResultTable <- tpptrSplineFitAndTest(data = fitData,
                                                  factorsH1 = "comparisonFactor",
                                                  resultPath = resultPath,
                                                  ggplotTheme = ggplotTheme,
                                                  doPlot = plotCurves,
                                                  splineDF = splineDF,# settings for spline fits
                                                  nCores = nCores, 
                                                  verbose = verbose,
                                                  additionalCols = annotData)
    if (flagDoWrite){
      save(list = c("splineFitResultTable"), 
           file = file.path(pathDataObj, "trResultsSplineFit.RData"))
    }
    resultTable <- left_join(resultTable, 
                             splineFitResultTable %>% 
                               mutate(Protein_ID = as.character(Protein_ID)),
                             by = "Protein_ID")
  }
  
  # Save result table to file if required
  if (flagDoWrite){
    save(list = c("resultTable"), 
         file = file.path(pathDataObj, "results_TPP_TR.RData"))
  }
  
  
  
  ## Save result table as xlsx spreadsheet:
  if (xlsxExport ){
    if (flagDoWrite){
      if (expNum > 1){
        ## Determine background colors for the columns belonging to the same experiment:
        compNums <- assignCompNumber_to_expName(compDF=expComps, expNames=expNames)
        expColors <- plotColors(expConditions = expConds, comparisonNums = compNums)
      } else {
        expColors <- NULL
      }      
      # Create output table:
      tppExport(tab=resultTable,  file=file.path(pathExcel, "results_TPP_TR.xlsx"), 
                expColors=expColors, expNames=expNames)
    } else {
      message("Cannot produce xlsx output because no result path is specified.")
    }
  }
  
  ## --------------------------------------------------------------------------------------------
  ## Create QC plots:
  if (flagDoWrite){
    pdf(file=file.path(resultPath, "QCplots.pdf"), width=8, height=9)
    
    ## 1. QC plot to illustrate median curve fits:
    message("Creating QC plots to visualize median curve fits...")
    if (normalize){
      qcPlotMedianFit <- normResults[["qcPlotObj"]]
      suppressWarnings(print(qcPlotMedianFit))
    }
    message("done.\n")
    
    ## 2. QC plot to correlate fold changes before and after normalization between all experiments:
    message("Creating QC plots to visualize normalization effects...")
    qcPlotCorrelateGroupsRaw  <- tppQCPlotsCorrelateExperiments(tppData=trData, 
                                                                annotStr="Non-normalized Fold Changes", 
                                                                ggplotTheme=ggplotTheme)
    if (normalize){
      qcPlotCorrelateGroupsNorm <- tppQCPlotsCorrelateExperiments(tppData=trDataNormalized, 
                                                                  annotStr="Normalized Fold Changes", 
                                                                  ggplotTheme=ggplotTheme)
    }
    for (pn in names(qcPlotCorrelateGroupsRaw)){
      suppressWarnings(print(qcPlotCorrelateGroupsRaw[[pn]]))
      if (normalize) suppressWarnings(print(qcPlotCorrelateGroupsNorm[[pn]]))
    }
    message("done.\n")
    
    if (any(methods == "meltcurvefit")){
      tpptrQCplots(resultTab = resultTable, expNames = expNames, 
                   expConditions = expConds, compDF = expComps, 
                   minR2 = pValFilter$minR2, ggplotTheme = ggplotTheme)
    }
    dev.off()
    ## --------------------------------------------------------------------------------------------
  }
  invisible(resultTable)
}
