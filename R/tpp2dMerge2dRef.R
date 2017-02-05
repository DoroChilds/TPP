#' @title Merge 2D-TPP result data with TPP-TR reference data
#'   
#' @description Merges 2D-TPP result data with TPP-TR reference data to generate
#' a big table including both results
#'  
#' @return A data frame with results merged from 2D-TPP and TPP-TR reference 
#' 
#' @param data DEPRECATED
#' @param trRef DEPRECATED
#' @param idVar DEPRECATED
#' @param resultTable_2D dataframe containing the 2D-TPP results
#' @param referenceDataSummary summarized reference data results. See details.
#' @param refIDVar character string indicating name of the columns containing 
#' the unique protein identifiers in the reference data set
#'   
#' @details \code{referenceSummary} contains summary statistics like median 
#' fold changes and is produced by the function
#'  \code{\link{tpp2dCreateTPPTRreference}}. It summarizes the results of a 
#'  TPP-TR analysis of a reference data set. 
#'  A reference data set is the a output of a TR experiment without drug 
#'  treatment on the same cell line as resultTable_2D. 
#' 
#' @examples  
#' data(panobinostat_2DTPP_smallExample)
#' config_tpp2d <- panobinostat_2DTPP_config
#' data_tpp2d <- panobinostat_2DTPP_data
#' tpp2dResults <- analyze2DTPP(configTable = config_tpp2d, 
#'                              data = data_tpp2d,
#'                              methods=c("doseResponse"),
#'                              createReport="none",
#'                              nCores=1,
#'                              idVar = "representative",
#'                              addCol = "clustername",
#'                              intensityStr = "sumionarea_protein_",
#'                              nonZeroCols = "qusm") 
#' 
#' trRef <- file.path(system.file("data", package="TPP"), 
#' "TPPTR_reference_results_HepG2.RData")
#' 
#' annotatedTable <- tpp2dMerge2dRef(resultTable_2D = tpp2dResults,
#'                                   referenceDataSummary = trRef)
#'
#' @seealso \code{\link{tpp2dCreateTPPTRreference}}
#' @export
#' 
tpp2dMerge2dRef <- function(resultTable_2D, 
                            referenceDataSummary,
                            refIDVar = "Protein_ID",
                            idVar = NULL,
                            data = NULL,
                            trRef = NULL){
  # Problem: is this function necessary? It is only invoked internally by 
  # analyze2DTPP to extract tppRefData$sumResTable$summary and merge it to the 
  # results of tpp2dCurveFit and tpp2dSplineFitAndTest
  # Additionally, there is also the function 'tpp2dAddAdditionalInfo' which 
  # is exported as well.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  tppRefData <- NULL
  
  if (!missing(data)){
    warning("`data` is deprecated.", call. = TRUE)
  }
  
  if (!missing(trRef)){
    warning("`trRef` is deprecated.", call. = TRUE)
  }
  
  if (!missing(idVar)){
    warning("`idVar` is deprecated.", call. = TRUE)
  }
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), 
                    c("resultTable_2D", "referenceDataSummary"))
  
  message("Function ", match.call()[[1]], ": Merging reference data summary to the 2D-TPP results.")
  
  # load tppRefData:
  if (is.data.frame(referenceDataSummary)){
    
    refDat <- referenceDataSummary
    
  } else if (is.character(referenceDataSummary)){
    
    if (file.exists(referenceDataSummary)){
      
      message("Function ", match.call()[[1]], ": Loading reference data from file ", referenceDataSummary)
      
      load(referenceDataSummary)
      refDat <- tppRefData$sumResTable$summary
      
    } else {
      
      stop("The path '", referenceDataSummary, "' that was given for the reference data summary does not exist.")
      
    }
    
  } else {
    
    stop("'referenceDataSummary' must be data frame created by the function 'tpp2dCreateTPPTRreference', or a path to such an object.")
  }
  
  # Obtain settings used for data import (stored as attribute of imported data):
  idVar <- checkAndReturnDataSetting(attr(resultTable_2D, "importSettings"), 
                                     "proteinIdCol", 
                                     colnames(resultTable_2D))
  
  if (is.null(attr(refDat, "importSettings"))){
    attr(refDat, "importSettings") <- list(proteinIdCol = refIDVar) # to do: remove later
  }
  refIDVar <- checkAndReturnDataSetting(dataSettings = attr(refDat, "importSettings"),
                                        "proteinIdCol",
                                        colnames(refDat))
  
  ## Remove all common column names except the ID vars:
  refCols <- colnames(refDat)
  dataCols <- colnames(resultTable_2D)
  inters <- intersect(refCols, dataCols) %>% setdiff(c(idVar, refIDVar))
  keepCols <- setdiff(refCols, inters)
  refDat2 <- refDat[, keepCols]
  
  ## Check if all ids are present in the reference data. Warn, if not. Stop if 
  ## no id matches (something probably went wrong in this case, e.g. wrong ID 
  ## column selected).
  datIDs <- resultTable_2D[[idVar]]
  refIDs <- refDat2[[refIDVar]]
  
  missingRefIDs <- setdiff(datIDs, refIDs)
  
  if (length(missingRefIDs) == length(datIDs)){
    
    txt1 <- "None of the identifies in the reference data summary matches those in the 2D-TPP results. Merging not possible."
    txt2 <- "Are the id columns assigned correctly?"
    txt3 <- paste0("Hint: Current id column for the 2D-TPP results = '", idVar, 
                   "'. Current id column in the reference data summary = '", 
                   refIDVar, "'.")
    stop(txt1, " ", txt2, "\n", txt3)
    
  } else if (length(missingRefIDs) > 0){
    
    txt1 <- "The following entries from the 2D-TPP results were not found in the reference data summary:"
    txt2 <- paste(missingRefIDs, collapse = ", ")
    txt3 <- "The corresponding rows have been filled with missing values."
    warning(txt1, "\n", txt2, ".\n", txt3, call. = TRUE)

  }
  
  # merge 2D with referenceDataSummary
  out <- merge(x = resultTable_2D, y = refDat2, by.x = idVar, 
               by.y = refIDVar, all.x = TRUE, all.y = FALSE)
  
  ## Avoid changes in ordering after merging:
  out2 <- out[match(out[[idVar]], resultTable_2D[[idVar]]),]
  
  # return resulting table
  return(out2)
  
}
