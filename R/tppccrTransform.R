#' @title Transform fold changes of TPP-CCR experiment
#' @description Transform fold changes of TPP-CCR experiment to prepare them for
#'   dose response curve fitting.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, data = hdacCCR_data)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' # Perform transformation:
#' tppccrTransformed <- tppccrTransform(data=tppccrNorm)
#' # Obtain transformed measurements per replicate:
#' transf_replicate1 <- tppccrTransformed$Panobinostat_1
#' transf_replicate2 <- tppccrTransformed$Panobinostat_2
#' # Inspect transformed data in replicate 1:
#' effects_replicate1 <- featureData(transf_replicate1)$compound_effect
#' newData_repl1 <- data.frame(exprs(transf_replicate1), 
#'                               Type=effects_replicate1)[!is.na(effects_replicate1),]
#'                               
#' @return List of expressionSet objects storing the transformed fold changes, 
#' as well as row and column metadata. In each expressionSet \code{S}, the fold changes
#'   can be accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. Isobaric labels and the corresponding concentrations are 
#'   returned by \code{S$label} and \code{S$concentration}.
#' 
#' @param data expressionSet object containing the data to be transformed.
#' @param fcCutoff cutoff for highest compound concentration fold change. 
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#' 
#' @details  
#' Only proteins with fold changes bigger than
#' \code{[fcCutoff * (1 - fcTolerance)} or smaller than 
#' \code{1/(fcCutoff * (1 - fcTolerance))]} will be used for curve fitting.
#' Additionally, the proteins fulfilling the fcCutoff criterion without 
#' tolerance will be marked in the output column \code{meets_FC_requirement}.
#' 
#' @export
tppccrTransform <- function(data, fcCutoff=1.5, fcTolerance=0.1) {
  dataListTransformed <- list()  
  for (expName in names(data)){
    message("Transforming dataset: ", expName)
    dTmp   <- data[[expName]]
    fcOrig <- exprs(dTmp)
    
    ## Mark proteins that are stabilized or destabilized by compound treatment
    fcMaxConc <- fcOrig[, which.max(dTmp$concentration)]
    
    ## 1. Use threshold without tolerance. Will be reported in result table.
    flagS_strict <- fcMaxConc >= fcCutoff
    flagD_strict <- fcMaxConc <= 1/fcCutoff
    featureData(dTmp)$meets_FC_requirement <- flagS_strict|flagD_strict
    
    ## 2. Repeat with tolerance to include proteins close to the threshold.
    flagS <- fcMaxConc >= fcCutoff * (1 - fcTolerance)
    flagD <- fcMaxConc <= 1/(fcCutoff * (1 - fcTolerance))
    
    cpdEffect        <- rep(NA_character_, nrow(dTmp))
    cpdEffect[flagS] <- "stabilized"
    cpdEffect[flagD] <- "destabilized"
    featureData(dTmp)$compound_effect <- cpdEffect

    fcNew <- matrix(NA_real_, nrow=nrow(dTmp), ncol=ncol(dTmp), 
                    dimnames=list(featureNames(dTmp), colnames(dTmp)))

    # Transform FCs of proteins stabilized by cpd treatment to 
    # fc=0 for DMSO and fc=1 for hightest cpd conc
    iS <- which(flagS)
    fcNew[iS,]   <- (fcOrig[iS  ,] - 1) / (fcMaxConc[iS] - 1)

    # Transform FCs of proteins destabilized by cpd treatment to 
    # fc=1 for DMSO and fc=0 for hightest cpd conc
    iD <- which(flagD)
    fcNew[iD,] <- (fcOrig[iD,] - fcMaxConc[iD])/ (1-fcMaxConc[iD])

    ## Store transformed FCs in exprs field:
    exprs(dTmp) <- fcNew
    
    ## Store transformed FCs in featureData so that it can be compared to the 
    ## untransformed values in the result table:
    fcNamesTransf <- paste(colnames(fcNew), "transformed", sep="_")
    pData(featureData(dTmp))[,fcNamesTransf] <- fcNew
    
    dataListTransformed[[expName]] <- dTmp
  }  
  message("Transformation complete.\n")
  return(dataListTransformed)  
}
