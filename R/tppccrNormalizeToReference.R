#' @title Normalize fold changes of TPP-CCR experiment to a reference column
#' @description Normalize fold changes of TPP-CCR experiment to a reference 
#' column (usually that with the lowest concentration) to ensure that the 
#' transformation by \link{tppccrTransform} yields values between 0 and 1.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, data = hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' # Normalize to lowest concentration (in the first column):
#' tppccrNormToRef <- tppccrNormalizeToReference(data=tppccrNorm, refCol=1)
#' # Obtain results per replicate:
#' refTransf_replicate1 <- tppccrNormToRef$Panobinostat_1
#' head(exprs(refTransf_replicate1))
#' # Perform transformation:
#' tppccrTransformed <- tppccrTransform(data=tppccrNormToRef)
#' # Obtain transformed measurements per replicate:
#' transf_replicate1 <- tppccrTransformed$Panobinostat_1
#' transf_replicate2 <- tppccrTransformed$Panobinostat_2
#' # Inspect transformed data in replicate 1:
#' effects_replicate1 <- featureData(transf_replicate1)$compound_effect
#' newData_repl1 <- data.frame(exprs(transf_replicate1), 
#'                               Type=effects_replicate1)[!is.na(effects_replicate1),]
#'                               
#' @return List of expressionSet objects storing the normalized fold changes, 
#' as well as row and column metadata. In each expressionSet \code{S}, the fold 
#' changes can be accessed by \code{exprs(S)}. Protein expNames can be accessed 
#' by \code{featureNames(S)}. Isobaric labels and the corresponding 
#' concentrations are returned by \code{S$label} and \code{S$concentration}.
#' 
#' @param data expressionSet object containing the data to be normalized
#' @param refCol column number to use as a reference. Will contain only 1s after 
#' the normalization.
#' 
#' @export
tppccrNormalizeToReference <- function(data, refCol=NULL){
  if (is.null(refCol)){
    refCol <- 1 # use lowest concentration by default.
  }
  dataListNormed <- list()  
  for (expName in names(data)){
    message("Normalizing dataset: ", expName, " to reference column ", refCol)
    dTmp   <- data[[expName]]
    foldChanges <- exprs(dTmp)
    refFC <- foldChanges[,refCol]
    refFC[which(refFC==0)] = 1e-15
    normedValues <- foldChanges / refFC
    exprs(dTmp) <- normedValues
    
    ## Store transformed FCs in featureData so that it can be compared to the 
    ## untransformed values in the result table:
    fcNamesRef <- paste(colnames(normedValues), "normalized_to_lowest_conc", sep="_")
    pData(featureData(dTmp))[,fcNamesRef] <- normedValues
    
    dataListNormed[[expName]] <- dTmp
  }
  return(dataListNormed)
}