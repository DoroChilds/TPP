#' @title Normalize data from TPP-CCR experiments
#' @description Normalize each fold change column by its median.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, data = hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' head(Biobase::exprs(tppccrNorm[[1]]))
#' 
#' @return List of expressionSet objects storing the normalized fold changes, as well as 
#'   row and column metadata. In each expressionSet \code{S}, the fold changes
#'   can be accessed by \code{Biobase::exprs(S)}. Protein names can be accessed by 
#'   \code{featureNames(S)}. Isobaric labels and the corresponding concentrations are 
#'   returned by \code{S$label} and \code{S$concentration}.
#'   
#' @param data list of expressionSets with measurements to be normalized
#' @export
tppccrNormalize <- function(data){
  if (!is.list(data) || is.data.frame(data)) {
    stop("'data' needs to be a list of expressionSets, in which the field names correspond to the experiment names.")
  }
  normData <- list()
  for (expName in names(data)){
    message("Normalizing dataset: ", expName)
    
    dTmp <- data[[expName]]
    
    ## Compute normalization coefficients:
    fcOld <- Biobase::exprs(dTmp)
    fcMedians <- apply(fcOld, 2, median, na.rm=TRUE)
    normCoeffs <- 1/fcMedians  
    
    ## Normalize using the computed coefficients:
    dTmpNorm <- applyCoeffs(data=dTmp, coeffs=normCoeffs)
    
    ## Store normalized fold changes in the featureData:
    fcNamesNorm <- paste(colnames(fcOld), "median_normalized", sep="_")
    pData(featureData(dTmpNorm))[,fcNamesNorm] <- Biobase::exprs(dTmpNorm)
    
    ## Store normalized dataset in output list:
    normData[[expName]] <- dTmpNorm
  }  
  message("Normalization complete.\n")
  return(normData)
}