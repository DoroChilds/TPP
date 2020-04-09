#' @title Create TPP-TR reference for 2D-TPP experiment
#'   
#' @description Performs a reference analysis of a TPP-TR experiment and generates boxplots for
#'  the distribution of fold changes at the different temperatures if desired.
#' 
#' @return A TPP-TR reference object for a certain cell line with different supporting files in a 
#'  desired output directory. The main object which is of interest for further analysis is the 
#'  \code{trRefData.RData} file. This is the file to which a referencing system path has to be 
#'  indicated when a function as \code{\link{tpp2dSplineFitAndTest}}
#'  require to input a TPP-TR reference object. 
#'  The RData file consists of list carrying four different items: \enumerate{ \item tppCfgTable: 
#'  the TPP-TR configtable which was used for generating this object \item sumResTable a list of
#'  two elements 1. detail: the exact result data from the TR analysis and 2. summary. a summary
#'  of the analyzed TR data comprising the median and standard deviation values of the measurements
#'  at the different temperatures (encoded by the isobaric labels) \item temperatures a table listing
#'  the temperatures which were used in the TR experiment in the different replicates \item lblsByTemp 
#'  a table matching each temperature to an isobaric label}
#' 
#' @param trConfigTable config file for a reference TR dataset
#' @param trDat list of dataframes, containing fold change measurements and 
#'   additional annotation columns to be imported. Can be used instead of 
#'   specifying the file path in the \code{configTable} argument.
#' @param resultPath character string containing a valid system path to which folder output files 
#'   will be written 
#' @param outputName character string which will be used as name of the output folder
#' @param createFCboxplots boolean flag indicating whether quality control boxplots are to be plotted
#' @param idVar character string indicating which column of the data table contains the unique
#'  protein ids
#' @param fcStr character string indicating which columns contain fold changes
#' @param qualColName character string indicating which column contain protein 
#'  identification quality measures 
#' @param normalize boolean argument stating whether the data should be normalized or not
#' 
#' @export
tpp2dCreateTPPTRreference <- function(trConfigTable=NULL, 
                                      trDat=NULL,
                                      resultPath=NULL, 
                                      outputName=NULL, 
                                      createFCboxplots=FALSE, 
                                      idVar="gene_name", 
                                      fcStr="rel_fc_", 
                                      qualColName="qupm", 
                                      normalize=TRUE){
  # set options
  options("TPPTR_plot" = FALSE)
  options("TPPTR_CI" = FALSE)
  
  if (is.null(trConfigTable)){
    stop("Please specify a valid argument for trConfigTable!")
  }
  
  if (is.null(resultPath)){
    stop("You have to specify a valid resultPath!")
  } else if (!file.exists(resultPath)){
    dir.create(resultPath, recursive = TRUE)
  }
  
  # create output name if outName is undefined 
  if(is.null(outputName) || !is.character(outputName)){
    stop("Please define a valid character string for the variable outputName!")
  }
  
  # generate outPath
  outPath = file.path(resultPath, outputName)
  if(file.exists(outPath)){
    outputName = paste(format(Sys.time(),'%Y-%m-%d_%H-%M'), "TPPTR_reference", 
                       sep="_")
    outPath = file.path(resultPath, outputName)
    message(paste("Changed outputName to ", outputName, 
                  ", to prevent unitentional overwriting of
                    previous results!"))
  }
  dir.create(outPath)
  message("Results will be written to ", outPath)
  
  # read in trConfigTable
  trConfigTable <- importFct_readConfigTable(trConfigTable)
  
  # determine temperatures
  temperatures <- trConfigTable[,which(!colnames(trConfigTable) %in% 
                                         c("Experiment", "Condition", "Path"))]
  lbls <- colnames(temperatures)
  temps <- temperatures[1, lbls]
  lblsByTemp <- data.frame(temp=as.character(temps), lbl=lbls)
  
  # define patterns of columns that are to be represented in the output table
  wantedColPatterns = c("meltPoint",
                        "slope",
                        "plateau",
                        "R_sq",
                        qualColName)
  
  # perform TPP-TR data import
  trData <- tpptrImport(configTable=trConfigTable, 
                        data=trDat,
                        idVar=idVar, 
                        fcStr=fcStr,
                        qualColName=qualColName)
  
  # normalize data if requested
  if (normalize){
    normResults <- tpptrNormalize(data=trData, 
                                  normReqs = tpptrDefaultNormReqs(), 
                                  qcPlotTheme = tppDefaultTheme(), 
                                  qcPlotPath = NULL, fixedReference = NULL)
    trDataNormalized <- normResults[["normData"]]
  }else {
    trDataNormalized <- trData
  }
  
  # write outout data file
  if (!is.null(outPath)){
    save(trDataNormalized, file=file.path(outPath, "normalizedData.RData"))    
  }
  
  # fit melting curves:
  trDataFitted <- tpptrCurveFit(data=trDataNormalized, dataCI=NULL, 
                                resultPath=outPath,
                                ggplotTheme=NULL, doPlot=FALSE,
                                startPars=c("Pl"=0, "a"=550, "b"=10), 
                                maxAttempts=500, 
                                nCores='max', verbose=FALSE)
  
  # analyze melting curves and create result table
  meltCurveResultTable <- tpptrAnalyzeMeltingCurves(data = trDataFitted) %>%
      rename(meltcurve_plot = plot) %>% 
      mutate(meltcurve_plot = as.character(meltcurve_plot)) %>%
      mutate(Protein_ID = as.character(Protein_ID))
  
  # save result table
  save(meltCurveResultTable, file=file.path(outPath, "resultTable.RData"))    
  
  # generate summary
  sumResTable = summarizeResultTable(meltCurveResultTable, wantedColPatterns, 
                                     temperatures, fcStr) 
  tppRefData = list(#'userCfgTable'=trConfigTable,
    'tppCfgTable'=trConfigTable,
    'temperatures'=temperatures,
    'sumResTable'=sumResTable,
    'lblsByTemp'=lblsByTemp,
    'idVar'=idVar,
    'fcStr'=fcStr,
    'qualColName'=qualColName)
  
  # save tppRefData
  save(tppRefData, file=file.path(outPath, "trRefData.RData"))
  
  # generate fold-change boxplots if desired
  if(createFCboxplots){
    # generate TPPTR reference object
    trRefObject <- tpp2dTRReferenceObject(tppRefPath=tppRefData)
    boxPlotPath = file.path(outPath, "fcBoxplots")
    dir.create(boxPlotPath)
    linkCol = rep("not_available", nrow(tppRefData$sumResTable$summary))
    # loop over all protein IDs and create boxplot
    plotList <- invisible(
      lapply(1:nrow(tppRefData$sumResTable$summary), function(row){
        protID = as.character(tppRefData$sumResTable$summary$Protein_ID[row])
        l_plot <- try(trRefObject$createFCBoxPlot(protID), silent=TRUE)
        # if no error occurred during plot generation, save plot
        if (class(l_plot) != "try-error") {
          fileName = paste("fcBoxpl_", gsub("\\.", "_", protID), ".pdf", sep="")
          savePath = file.path(boxPlotPath, fileName)
          message(paste("Saving ", fileName, "...", sep=""))
          try(ggsave(filename=savePath, plot=l_plot, width=9, height=9))
          linkCol[row] = file.path(basename(boxPlotPath), fileName)
          return(l_plot)
        } else{
          return(NULL)
        }
      }))
    
    tppRefData$sumResTable$detail["fc_Boxplot"] = linkCol
    tppRefData$sumResTable$summary["fc_Boxplot"] = linkCol
  }
  
  # write output files
  writeTRRefOutputTable(list('Configuration'=trConfigTable, 
                             'Summary'=tppRefData$sumResTable$summary, 
                             'Details'=tppRefData$sumResTable$detail), 
                        outPath, outputName, type='xlsx', 
                        multiSheet=TRUE)
  writeTRRefOutputTable(list('Configuration'=trConfigTable, 
                             'Summary'=tppRefData$sumResTable$summary, 
                             'Details'=tppRefData$sumResTable$detail), 
                        outPath, outputName, type='txt')
}
