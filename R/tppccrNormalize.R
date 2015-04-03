#' @title Normalize data from TPP-CCR experiments
#' @description Normalize each fold change column by its median.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config_repl1, data = hdacCCR_data_repl1)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' head(exprs(tppccrNorm))
#' 
#' @return ExpressionSet object storing the normalized fold changes, as well as 
#'   row and column metadata. In each ExpressionSet \code{S}, the fold changes
#'   can be accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. TMT labels and the corresponding temperatures are 
#'   returned by \code{S$labels} and \code{S$temperatures}.
#'   
#' @param data expressionSet with measurements to be normalized
#' @export
tppccrNormalize <- function(data){
  message("Normalizing data ...")
  
  ## Compute normalization coefficients:
  fcOld <- exprs(data)
  fcMedians <- apply(fcOld, 2, median, na.rm=TRUE)
  normCoeffs <- 1/fcMedians  
  
  ## Normalize using the computed coefficients:
  dataNew <- applyCoeffs(data=data, coeffs=normCoeffs)
  
  ## Return result
  message("done.")
  return(dataNew)
}