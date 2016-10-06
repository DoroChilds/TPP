#' @title Analyze fitted curve parameters to detect significant shifts in 
#'   melting points.
#' @description Compute p-values for the pairwise comparisons of melting curve 
#'   shifts between different conditions.
#' @details The \code{pValParams} argument is a list that can contain optional parameters 
#'   for the chosen p-value computation \code{pValMethod}. The following options are
#'   available: \enumerate{ \item \code{pValMethod = "robustZ"}:
#'   
#'   \code{pValParams=list(binWidth=[your_binWidth])}.}
#' @return A data frame in which the fit results are stored row-wise for each 
#'   protein.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(hdacTR_config, hdacTR_data)
#' tpptrNorm <- tpptrNormalize(data=tpptrData, 
#'                             normReqs=tpptrDefaultNormReqs())
#' normalizedData <- tpptrNorm$normData
#' \dontrun{
#' # Fit melting curves to each protein 
#' # (can take some time depending on device used):
#'  fittedData <- tpptrCurveFit(normalizedData, nCores=1)
#'  resultTable <- tpptrAnalyzeMeltingCurves(fittedData)
#'  subset(resultTable, fulfills_all_4_requirements)$Protein_ID
#' }
#' 
#' @param data list of ExpressionSets containing fold changes and metadata. Their
#'   featureData fields contain the fitted melting curve parameters.
#' @param pValMethod Method for p-value computation. Currently restricted to
#'   'robustZ' (see Cox & Mann (2008)).
#' @param pValParams optional list of parameters for p-value computation.
#' @param pValFilter optional list of filtering criteria to be applied before 
#'   p-value computation.
#' @references Cox, J., & Mann, M. (2008). MaxQuant enables high peptide 
#'   identification rates, individualized ppb-range mass accuracies and 
#'   proteome-wide protein quantification. Nature biotechnology, 26(12), 
#'   1367-1372.
#' @export
tpptrAnalyzeMeltingCurves <- function(data, pValMethod="robustZ", 
                                      pValFilter=list(minR2=0.8, maxPlateau=0.3),
                                      pValParams=list(binWidth=300)){
  message("Starting melting curve analysis.")
  expInfo <- sapply(data, annotation)
  
  dataSplit   <- retrieveDataFromESets_TR(data)
  curveParDF  <- dataSplit$curveParDF
  compDF      <- createComparisonTable(infoTable=expInfo)
  
  if (!is.null(compDF)){
    pValCols      <- pValFctPerformComparisons(curveParsAllExp = curveParDF, 
                                               method = pValMethod, 
                                               controlFilter = pValFilter, 
                                               controlpVal = pValParams, 
                                               comparisons = compDF)
    qualCheckCols <- checkResultCols_tpptr(pValDF=pValCols, 
                                           curveParDF=curveParDF, 
                                           comparisons=compDF)    
  } else{
    pValCols <- qualCheckCols <- NULL
  }
  resultTable <- mergeOutputTables_TR(dataList=dataSplit, pValDF=pValCols, 
                                      qualCheckDF=qualCheckCols)
  
  return(resultTable)
}