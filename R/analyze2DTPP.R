#' @title Analyze 2D-TPP experiment
#'   
#' @description Performs analysis of a 2D-TPP experiment by invoking routines 
#'   for data import, data processing, fold change computation, median normalization, 
#'   TPP-CCR curve fitting, plotting and production of the result table.
#'  
#' @return A data frame in which the fit results are stored row-wise for each
#'   protein at the different temperatures.
#'   
#' @references Becher, I., Werner, T., Doce, C., Zaal, E. A., Berkers, C. R., T"ogel, I., 
#' Salzer, E., Bantscheff, M., Savitski, M. M. (2016) Comprehensive thermal and chemoproteomics 
#' profiling identifies phenylalanine hydroxylase as a potent off-target of the histone 
#' deacetylase inhibitor panobinostat. Nature Chemical Biology (accepted)
#' 
#' @details Invokes the following steps: \enumerate{ \item Import data using the
#'   \code{\link{tpp2dImport}} function. \item Remove zero sumionarea values. 
#'   \item Compute fold changes from raw data (sumionarea) 
#'   \item Perform normalization by fold 
#'   change medians (optional) using the \code{\link{tpp2dNormalize}} function.
#'   To perform normalization, set argument \code{normalize=TRUE}.}
#'
#' @examples 
#' data("panobinostat_2DTPP_smallExample")
#' config_tpp2d <- panobinostat_2DTPP_config
#' data_tpp2d <- panobinostat_2DTPP_data
#' tpp2dResults <- analyze2DTPP(configFile = config_tpp2d, 
#'                              data = data_tpp2d,
#'                              fcStr = NULL,
#'                              methods=c("doseResponse"),
#'                              createReport="none",
#'                              nCores=1)
#' 
#' @param configFile dataframe, or character object with the path to a file, 
#'   that specifies important details of the 2D-TPP experiment. See Section 
#'   \code{details} for instructions how to create this object.
#' @param data single dataframe, containing fold change measurements and 
#'   additional annotation columns to be imported. Can be used instead of 
#'   specifying the file path in the \code{configTable} argument.
#' @param resultPath location where to store dose-response curve plots and 
#'   results table.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the prefix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the prefix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param methods vector of character strings that indicate which methods should be used 
#'   for the analysis (default: c("doseRespone"), alternative: c("splineFit") or 
#'   c("doseRespone", "splineFit")) 
#' @param addCol vector of chracter strings indicating which additional columns to include 
#'   from the input data 
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param compFc boolean flag which indicates whether to perform fold change computation 
#'   regarding reference column from sumionareas (default: TRUE)
#' @param normalize perform median normalization (default: TRUE).
#' @param nCores either a numerical value given the desired number of CPUs, or 
#'   'max' to automatically assign the maximum possible number (default).
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param xlsxExport produce results table in xlsx format and store at the 
#' location specified by the \code{resultPath} argument.
#' @param plotAll boolan value indicating whether all dose response curves should
#' be generated. Deactivating plotting decreases runtime.
#' @param plotAllR2 boolan value indicating whether all dose response curves which
#' fullfill the demanded criterias (Rsquared, maximum plateau) should be generated. 
#' Deactivating plotting decreases runtime.
#' @param plotSingle boolan value indicating whether all dose response curves which
#' fullfill the demanded criterias (Rsquared, maximum plateau) should be generated. 
#' Deactivating plotting decreases runtime.
#' @param fractAbund boolean variable, if set to TRUE additional information concerning
#' sumionarea fractional abundance and dmso1 vs. dmso2 of adjacent temperatures is 
#' added to the output table
#' @param addInfo boolean variable, if set to TRUE additional information on counts of
#' stabilization and destabiliazation of each protein is added to the output table
#' @param trRef character string containing a valid system path to a previously generated TPP-TR
#' reference object
#' @param refFcStr character string indicating which columns in the reference data set contain 
#' the fold change values
#' @param createReport character string indicating whether a markdown report should be created
#'  and which format it have (default: "html_document", alternative: "pdf_document" or "none")
#' 
#' @export
analyze2DTPP <- function(configFile = NULL, data = NULL, 
                         resultPath = NULL, idVar = "representative", 
                         fcStr = "rel_fc_protein_", 
                         intensityStr = "sumionarea_protein_",   
                         naStrs = c("NA", "n/d", "NaN", "<NA>"), 
                         methods = c("doseRespone", "splineFit"),
                         qualColName="qupm", compFc=TRUE, normalize=TRUE, addCol=NULL,
                         nCores="max", nonZeroCols="qupm", fcTolerance=0.1, r2Cutoff=0.8,  
                         fcCutoff=1.5, slopeBounds=c(1,50), fractAbund=FALSE, xlsxExport=FALSE, 
                         plotAll=FALSE,plotAllR2=FALSE, plotSingle=FALSE, trRef=NULL, 
                         refFcStr="norm_rel_fc_protein_", addInfo=FALSE, createReport="html_document") {
  
  message("This is TPP version ", packageVersion("TPP"),".")
  
  # check configTable and read in table if necessary
  configTable <- tpp2dEvalConfigTable(configFile)
  
  # import data
  Data2d <- tpp2dImport(configTable=configTable, data=data, 
                            idVar=idVar, addCol=addCol, intensityStr=intensityStr, 
                            qualColName=qualColName, fcStr=fcStr)
  
  # compute fold changes if requested
  if (compFc){
    if (is.null(fcStr)){
      fcStr <- "rel_fc_protein_"
    }
    Data2d <- tpp2dComputeFoldChanges(configTable=configTable, data=Data2d, 
                                      intensityStr=intensityStr) 
    
  }else if (is.null(fcStr)){
    stop("Fold changes need to either be supplied in the input data by specifying a prefix in fcStr 
         or they need to be computed be setting comFc to TRUE!")
  }
  
  # do median normalization of fold changes 
  if (normalize){
    NormData2d <- tpp2dNormalize(configTable=configTable, data=Data2d, fcStr=fcStr)
    
    # Make sure the TPP-CCR routine uses the correct columns, when there was 
    # normalization before:
    fcStrUpdated <- paste("norm", fcStr, sep="_")
  }else{
    NormData2d <- Data2d
    fcStrUpdated <- fcStr
  }
  
  # filter out row with no quality information
  if (length(which(is.na(NormData2d[[qualColName[1]]])))!=0){
    NormData2d <- NormData2d[-which(is.na(NormData2d[[qualColName[1]]])),]
  }
  
  # filter out proteins that have duplicated unique_ID to prevent crash 
  if (length(which(duplicated(NormData2d$unique_ID)))!=0){
    message("There are duplicated proteins in your experimental conditions! These are filtered out to run this analysis!
Please check your data quality and consider pre-filtering!")
    NormData2d <- NormData2d[-which(duplicated(NormData2d$unique_ID)),]
  }
  
  # calculate fractioncal abundance per curve
  if (fractAbund){
    NormData2d <- tpp2dCalcFractAbundance(configTable=configTable, data=NormData2d, 
                                          intensityStr=intensityStr, idVar=idVar)
  }
  
  if ("doseResponse" %in% methods){
    # create CCR config file list
    CCR2dConfig <- tpp2dCreateCCRConfigFile(configTable=configTable)
    
    # run TPP-CCR
    analysisResults <- tpp2dCurveFit(configFile = CCR2dConfig, 
                                      data = NormData2d, 
                                      nCores = nCores, 
                                      fcStr = fcStrUpdated, 
                                      idVar = "unique_ID", 
                                      nonZeroCols = nonZeroCols, 
                                      r2Cutoff = r2Cutoff, 
                                      fcCutoff = fcCutoff, 
                                      slopeBounds = slopeBounds,
                                      fcTolerance = fcTolerance)
    
    if (plotAll){
      # generate joint plots for all proteins detected
      plotList <- tpp2dPlotCCRAllCurves(configTable=configTable, data=analysisResults, 
                                        idVar=idVar, fcStr=fcStrUpdated)
      # write output file with plots
      tpp2dExportPlots(plotList=plotList, resultPath=resultPath, type="all")
    }
    if (plotAllR2){
      # generate joint plots for all proteins detected with sufficient R2
      plotGoodList <- tpp2dPlotCCRGoodCurves(configTable=configTable, data=analysisResults, 
                                             idVar=idVar, fcStr=fcStrUpdated)
      # write output file with plots
      tpp2dExportPlots(plotList=plotGoodList, resultPath=resultPath, type="good")
    }
    if (plotSingle){
      # generate single plots for all protein in each condition fitted with sufficient R2
      plotSingleList <- tpp2dPlotCCRSingleCurves(configTable=configTable, data=analysisResults, 
                                                 idVar=idVar, fcStr=fcStrUpdated)
      # write output file with plots
      tpp2dExportPlots(plotList=plotSingleList, resultPath=resultPath, type="single")
    }
  }else{
    analysisResults <- NormData2d
  }
  
  # do spline fit over tpp-tr reference
  if (("splineFit" %in% methods) && !is.null(trRef)){
    # do f-test for splines fit
    analysisResults <- tpp2dSplineFitAndTest(data_2D = analysisResults, 
                                             trRefDataPath = trRef, 
                                             idVar = idVar,
                                             fcStr = fcStrUpdated, 
                                             refFcStr = refFcStr,
                                             resultPath = resultPath,
                                             verbose = verbose)
    
  }else if("splineFit" %in% methods){
    message("The spline fit and corresponding f-Test could not be performed, as no TPP-TR reference dataset was specified!
Please check the file path you have specified for trRef!")
  }
  
  
  
  # add additional information e.g. how often protein was stabilized/destabilized if desired
  if (addInfo){
    analysisResults <- tpp2dAddAdditionalInfo(data = analysisResults)
  }
  
  # add TR reference columns to result table
  if (!is.null(trRef) && file.exists(trRef)){
    analysisResults <-tpp2dMerge2dRef(data = analysisResults, 
                                      trRef = trRef, idVar = idVar)
  }
  
  # export results
  if (!is.null(resultPath)){
    tpp2dExport(configTable = configTable, tab = analysisResults, 
                resultPath = resultPath, 
                idVar = "Protein_ID", fcStr = fcStr, intensityStr = intensityStr, # new: idVar = "Protein_ID"
                addCol = addCol, 
                normalizedData = normalize, trRef = trRef) 
  }
  
  # create markdown report
  if ((createReport!="none") && !is.null(resultPath)){
    tpp2dCreateReport(resultPath = resultPath, 
                      configFile = configFile, 
                      normalize = normalize,
                      configTable = configTable, 
                      resultTable = analysisResults, 
                      idVar = "Protein_ID", 
                      fcStr = fcStr, 
                      fcStrUpdated = fcStrUpdated, 
                      documentType = createReport,
                      intensityStr = intensityStr, 
                      addCol = addCol, 
                      fcTolerance = fcTolerance, 
                      r2Cutoff = r2Cutoff, 
                      fcCutoff = fcCutoff, 
                      slopeBounds = slopeBounds, 
                      trRef = trRef)
  }
  
  # return output table
  return(analysisResults)
}
