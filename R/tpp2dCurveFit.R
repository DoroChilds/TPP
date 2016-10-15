#' @title Run TPP-CCR analysis for 2D-TPP experiment
#'   
#' @description Performs analysis of a TPP-CCR experiment by invoking the routine 
#'   for TPP-CCR curve fitting for each temperature of the sample.
#'  
#' @return A data frames in which the fit results are stored row-wise for each
#'   protein.
#'   
#' @examples 
#'   load(system.file("example_data/2D_example_data/referenceCCRConfig.RData", package="TPP"))
#'   load(system.file("example_data/2D_example_data/exampleRunCCRInput.RData", package="TPP"))
#'   CCRresults <- tpp2dCurveFit(configFile=exampleCCRConfig, data=exampleRunCCRInput, 
#'                                idVar="unique_ID")
#'   
#' @param configFile list of dataframes, that specifies important details of the 2D-TPP 
#'   experiment for each temperature. 
#' @param data data.frame, that contains the data of the 2D-TPP 
#'   experiment for each temperature. 
#' @param nCores numeric value stating how many cores are to be used for computation
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#'   
#' @export
tpp2dCurveFit <- function(configFile, data, nCores=1, 
                           naStrs=c("NA", "n/d", "NaN", "<NA>"), 
                           fcStr="norm_rel_fc_protein_", 
                           idVar=NULL, nonZeroCols="qupm",
                           r2Cutoff=0.8,  fcCutoff=1.5, slopeBounds=c(1,50),
                           fcTolerance=0.1){
  
  if (is.null(configFile) || is.null(data)){
    stop("Please specifiy valid dataframes for the arguments configTable and data!")
  }else if (is.null(idVar) || !is.character(idVar)){
    stop("Please specify a valid character string for idVar!")
  }else if (!idVar %in% colnames(data)){
    stop("Please specify an idVar character string argument that represents a suffix of one of 
         the column names of your data!")
  }else if (length(data[[idVar]])!=length(unique(data[[idVar]]))){
    stop("Please indicate an idVar character string that matches a column with unique identifiers!")
  }else{
    message(paste("Performing TPP-CCR dose response curve fitting and generating result table...", 
                  sep=" "))
    CCRresult <- suppressMessages(analyzeTPPCCR(configTable=configFile, 
                                                data=as.data.frame(data), nCores=nCores, 
                                                resultPath=NULL, plotCurves=FALSE, fcStr=fcStr, 
                                                naStrs=naStrs, xlsxExport=FALSE,idVar=idVar, 
                                                nonZeroCols=nonZeroCols, normalize=FALSE, 
                                                r2Cutoff=r2Cutoff, fcCutoff=fcCutoff, 
                                                slopeBounds=slopeBounds,fcTolerance=fcTolerance))
    # reformat colnames
    compound <- as.character(configFile$Experiment)
    colnames(CCRresult) <- sub(paste("*_", compound, sep=""), "", colnames(CCRresult))
    message("Done.")
    return(CCRresult) 
  }
}
