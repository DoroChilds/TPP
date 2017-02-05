#' @title Compute 2D-TPP fold changes
#' @description Computes fold changes by calculating fold changes of the sumionarea 
#'  relative to the reference column.
#'
#' @return A data.frame with additional columns with constitute fold changes calculated with 
#'  respect to the intensity values of the zero treatment column 
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
#' # View attributes of imported data (experiment infos and import arguments):
#' attr(datIn, "importSettings") %>% unlist
#' attr(datIn, "configTable")
#' 
#' # Compute fold changes:
#' datFC <- tpp2dComputeFoldChanges(data = datIn)
#' 
#' # View updated attributes. Now contain field 'fcStrNorm' indicating prefix
#' # of the fold change columns after normalization.
#' attr(datFC, "importSettings")["fcStr"]
#'                                  
#' @param data dataframe that contain the data for the 2D-TPP experiment
#' @param configTable DEPRECATED
#' @param intensityStr DEPRECATED
#' @param fcStr DEPRECATED
#' @param newFcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{newFcStr} will be pasted in front of
#'   the names of the experiments.
#' 
#'   
#' @export
tpp2dComputeFoldChanges <- function(configTable = NULL, 
                                    data, 
                                    intensityStr = NULL, 
                                    fcStr = NULL,
                                    newFcStr = "rel_fc_"){

  if (!missing(configTable)){
    warning("`configTable` is deprecated.", call. = TRUE)
  }
  
  if (!missing(intensityStr)){
    warning("`intensityStr` is deprecated.", call. = TRUE)
  }
  
  if (!missing(fcStr)){
    warning("`fcStr` is deprecated.", call. = TRUE)
  }
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data"))
  
  # Obtain config table used for data import (stored as attribute of imported data):
  configTable <- attr(data, "configTable")
  
  # Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(data, "importSettings")
  
  intensityStr <- importSettings$intensityStr
  
  # Check whether intensitStr has class character.
  if (is.null(intensityStr)){
    stop("attr(data, 'importSettings') must contain a field named 'intensityStr'.")
  } else if (!is.character(intensityStr)){
    stop("attr(data, 'importSettings')$intensityStr must be of class character.")
  } else {
    message("Looking for intensity column prefix: '", intensityStr, "'")
  }
  
  # Determine reference colnames for experiments
  intensity.col.ids <- grep(intensityStr, colnames(data))
  
  if (length(intensity.col.ids) == 0){
    stop("The given prefix for intensity columns (specified by attr(data, 'importSettings')$intensityStr) was not found in the column names of 'data'.")
  }
  
  message("Computing fold changes...")
  
  fc.cols <- sapply(intensity.col.ids, function(sc){
    return(paste(newFcStr, sub(intensityStr, "", colnames(data)[sc]), sep=""))
  })
  
  ref.col.id <- grep(paste(intensityStr, "0$", sep=""), colnames(data))
  if (length(ref.col.id)==0){
    ref.col.id <- grep(paste(intensityStr, "0.0$", sep=""), colnames(data))
    if (length(ref.col.id)==0){
      ref.col.id <- grep(paste(intensityStr, "0.00$", sep=""), colnames(data))
      if (length(ref.col.id)==0){
        ref.col.id <- grep(paste(intensityStr, "0.000$", sep=""), colnames(data))
      }
    }
  }
  
  ref.vals <- as.numeric(as.character(data[,ref.col.id])) 
  
  # Calculate fold change of respective values by dividing by the reference value
  fc.m <- data.frame(matrix(as.numeric(as.character(unlist(data[,intensity.col.ids]))),
                            nrow=nrow(data[,intensity.col.ids])) / ref.vals)
  
  names(fc.m) <- fc.cols
  data[fc.cols] <- fc.m
  
  # Update attributes: Store column prefix of the new fold change columns
  importSettings$fcStr <- newFcStr
  
  attr(data, "importSettings") <- importSettings
  
  message("Done.")
  
  return(data)
  
} 
