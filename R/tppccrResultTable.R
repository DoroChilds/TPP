#' @title Summarize results of a TPP-CCR study
#' @description \code{tppccrResultTable} summarizes the 
#' outcomes of a TPP-CCR study in a results table and includes quality information
#'   about the estimated dose response curves.
#' @param data list of expressionSet objects containing protein fold changes, as 
#' well as fitted curve parameters.
#' @param r2Cutoff quality criterion on dose response curve fit.
#' 
#'  @details \code{data} is a list of expressionSet objects created by
#' \code{\link{tppccrCurveFit}} or \code{\link{tppccrPlotCurves}}. 
#' It contains the isobaric labels and administered drug concentrations in the 
#' \code{phenoData} and user-defined protein properties (including dose response 
#' curve parameters) in the \code{featureData}. 
#' Protein IDs are stored in the \code{featureNames}.
#' 
#' If \code{data} is the output of \code{\link{tppccrPlotCurves}}, 
#' plot locations are given in the \code{plot} column of the \code{featureData}.
#' 
#'  
#' @return A data frame in which the results are stored row-wise for each
#'   protein, together with the original annotation from the input files.
#'
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, 
#'                            data=hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' tppccrTransformed <- tppccrTransform(data=tppccrNorm)
#' tppccrFitted <- tppccrCurveFit(data=tppccrTransformed, nCores=1)
#' tppccrResults <- tppccrResultTable(data=tppccrFitted)
#' subset(tppccrResults, passed_filter_Panobinostat_1 & passed_filter_Panobinostat_2)
#' 
#' @seealso \code{\link{tppccrCurveFit}},\code{\link{tppccrPlotCurves}}
#'  
#' @export
tppccrResultTable <- function(data, r2Cutoff=0.8){
  dataSplit <- retrieveDataFromESets_CCR(data = data)
  
  curveParDF  <- dataSplit$modelPars
  fcFilterDF  <- dataSplit$transfDF
  qualCheckCol <- checkResultCols_tppccr(curveParDF=curveParDF, 
                                         fcFilterDF=fcFilterDF, minR2=r2Cutoff)
  dataSplit$transfDF <- cbind(dataSplit$transfDF, qualCheckCol)
  
  outTable <- mergeOutputTables_CCR(dataList=dataSplit, qualCheckDF=qualCheckCol)
  
  return(outTable)  
}

