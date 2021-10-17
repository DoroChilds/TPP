#' @title Run TPP-CCR analysis for 2D-TPP experiment
#'   
#' @description Performs analysis of a TPP-CCR experiment by invoking the routine 
#'   for TPP-CCR curve fitting for each temperature of the sample.
#'  
#' @return A data frames in which the fit results are stored row-wise for each
#'   protein.
#'   
#' @examples 
#' # Preparation:
#' data(panobinostat_2DTPP_smallExample)
#' 
#' # Import data:
#' datIn <- tpp2dImport(configTable = panobinostat_2DTPP_config,
#'                       data = panobinostat_2DTPP_data,
#'                       idVar = "representative",
#'                       addCol = "clustername",
#'                       intensityStr = "sumionarea_protein_",
#'                       nonZeroCols = "qusm")
#' 
#' # Compute fold changes:
#' datFC <- tpp2dComputeFoldChanges(data = datIn)
#'
#' # Perform median normalization:
#' datNorm <- tpp2dNormalize(data = datFC)
#' 
#' # View updated attributes. Now contain field 'fcStrNorm' indicating prefix
#' # of the fold change columns after normalization.
#' attr(datNorm, "importSettings")["fcStrNorm"]
#' 
#' # Perform dose response curve fitting and pEC50 calculation:
#' datFit <- tpp2dCurveFit(data = datNorm)
#'   
#' @param configFile DEPCRECATED
#' @param data data frame that contains the data of the 2D-TPP 
#'   experiment for each temperature. 
#' @param nCores numeric value stating how many cores are to be used for computation
#' @param naStrs DEPCRECATED
#' @param fcStr DEPCRECATED
#' @param idVar DEPCRECATED
#' @param nonZeroCols DEPCRECATED
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change.
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting.
#' @param fcTolerance tolerance for the fcCutoff parameter. See details.
#'   
#' @export
tpp2dCurveFit <- function(configFile = NULL, 
                          data, 
                          nCores = 1, 
                          naStrs = NULL, 
                          fcStr = NULL, 
                          idVar = NULL, 
                          nonZeroCols = NULL,
                          r2Cutoff = 0.8,  
                          fcCutoff = 1.5, 
                          slopeBounds = c(1,50),
                          fcTolerance = 0.1){
  
  if (!missing(configFile)){
    warning("`configFile` is deprecated.", call. = TRUE)
  }
  
  if (!missing(naStrs)){
    warning("`naStrs` is deprecated.", call. = TRUE)
  }
  
  if (!missing(fcStr)){
    warning("`fcStr` is deprecated.", call. = TRUE)
  }
  
  if (!missing(idVar)){
    warning("`idVar` is deprecated.", call. = TRUE)
  }
  
  if (!missing(nonZeroCols)){
    warning("`nonZeroCols` is deprecated.", call. = TRUE)
  }
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data"))
  
  
  # Obtain config table used for data import (stored as attribute of imported data):
  configTable <- attr(data, "configTable")
  
  # Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(data, "importSettings")
  
  uniqueIdCol <- importSettings$uniqueIdCol
  
  # Check whether uniqueIdCol has class character.
  if (is.null(uniqueIdCol)){
    stop("attr(data, 'uniqueIdCol') must contain a field named 'uniqueIdCol'.")
  } else if (!is.character(uniqueIdCol)){
    stop("attr(data, 'importSettings')$uniqueIdCol must be of class character.")
  } else {
    message("Looking for unique ID column: '", uniqueIdCol, "'")
  }
  
  if (!uniqueIdCol %in% colnames(data)){
    stop("Please specify an uniqueIdCol character string argument that represents a suffix of one of 
         the column names of your data!")
  } else if (length(data[[uniqueIdCol]])!=length(unique(data[[uniqueIdCol]]))){
    stop("Please indicate an uniqueIdCol character string that matches a column with unique identifiers!")
  }
  
  nonZeroCols <- importSettings$nonZeroCols
  
  # Check whether nonZeroCols are valid column names.
  if (is.null(nonZeroCols)){
    stop("attr(data, 'importSettings') must contain a field named 'nonZeroCols'.")
  } else if (!is.character(nonZeroCols)){
    stop("attr(data, 'importSettings')$nonZeroCols must be of class character.")
  } else {
    message("Looking for nonZeroCols: '", nonZeroCols, "'")
  }
  
  if (!all(nonZeroCols %in% colnames(data))){
    stop("The given QC columns (specified by attr(data, 'importSettings')$nonZeroCols) were not found in the column names of 'data'.")
  }
  
  
  # Choose correct fold change column prefix (automatically detects whether
  # to use the prefix for normalized columns).
  finalFcPrefix <- obtain_fcStr_from_df_annotation(dat = data)
  
  message("Performing TPP-CCR dose response curve fitting and generating result table...")
  
  cfgIn <- convert_2D_cfgTable_to_CCR_cfgTable(configTable = configTable)
  
  datIn <- as.data.frame(data)
  
  CCRresult <- suppressMessages(
    analyzeTPPCCR(configTable = cfgIn, 
                  data = datIn, 
                  nCores = nCores, 
                  resultPath = NULL, 
                  plotCurves = FALSE, 
                  fcStr = finalFcPrefix, 
                  naStrs=c("NA", "n/d", "NaN", "<NA>"),
                  qualColName="qupm",
                  xlsxExport = FALSE, 
                  idVar = uniqueIdCol, 
                  nonZeroCols = nonZeroCols, 
                  normalize = FALSE, 
                  r2Cutoff = r2Cutoff, 
                  fcCutoff = fcCutoff, 
                  slopeBounds = slopeBounds,
                  fcTolerance = fcTolerance,
                  ggplotTheme = NULL,
                  verbose=FALSE)
  )
  
  # Remove compound name suffix from each column of TPPCCR output
  compound <- as.character(cfgIn$Experiment)
  allCols <- colnames(CCRresult)
  newCols <- sub(paste("*_", compound, sep=""), "", allCols)
  colnames(CCRresult) <- newCols
  message("Done.")
  
  # Transfer attributes to newly created data frame
  importSettings$r2Cutoff <- r2Cutoff  
  importSettings$fcCutoff <- fcCutoff
  importSettings$slopeBounds <- slopeBounds
  importSettings$fcTolerance <- fcTolerance
  importSettings$uniqueIdCol <- "Protein_ID"
  
  attr(CCRresult, "importSettings") <- importSettings
  attr(CCRresult, "configTable")    <- configTable
  
  return(CCRresult) 
}
