#' @title Median normalization of protein fold changes of 2D-TPP data
#'   
#' @description Normalizes fold changes retrieved from 2D-TPP experiment by dividing by the median fold
#'   change 
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
#' @param data data frame that contains the data for the 2D-TPP experiment
#' @param configTable DEPRECATED
#' @param fcStr DEPRECATED
#'   
#' @return A dataframe identical to the input dataframe except that the columns containing the
#'   fold change values have been normalized by their median.
#'   
#' @export 
tpp2dNormalize <- function(configTable = NULL, data, fcStr = NULL){
  
  ## ----------------------------------------------------------------------- ##
  ## Preparation
  ## ----------------------------------------------------------------------- ##
  
  if (!missing(configTable)){
    warning("`configTable` is deprecated.", call. = TRUE)
  }
  
  if (!missing(fcStr)){
    warning("`fcStr` is deprecated.", call. = TRUE)
  }
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  median_per_temp_and_conc = y = yNormalized = columnName <- NULL
  
  ## Check for missing function arguments
  checkFunctionArgs(match.call(), c("data"))
  
  ## Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(data, "importSettings")
  idVar <- checkAndReturnDataSetting(importSettings, "proteinIdCol", colnames(data))
  fcStr <- checkAndReturnDataSetting(importSettings, "fcStr", colnames(data))
  
  ## Define old and new measurement columns
  fcStrNorm <- paste("norm", fcStr, sep="_")
  fcCols <- grep(fcStr, colnames(data), value = TRUE)
  fcNormCols <- gsub(fcStr, fcStrNorm, fcCols)
  
  ## ----------------------------------------------------------------------- ##
  ## Perform normalization
  ## ----------------------------------------------------------------------- ##
  message("Performing median normalization per temperature...")
  
  dataLong <- data %>% subset(select = c(idVar, "temperature", fcCols)) %>%
    gather_("columnName", "y", gather_cols = fcCols) %>%
    mutate(y = as.numeric(as.character(y)))
  
  normCoeffs <- dataLong %>% 
    group_by_(.dots = c("temperature", "columnName")) %>%
    summarise(median_per_temp_and_conc = median(y, na.rm = TRUE))
  
  dataNormed <- dataLong %>% 
    left_join(normCoeffs, by = c("temperature", "columnName")) %>% 
    mutate(yNormalized = y / median_per_temp_and_conc) %>%
    select(-median_per_temp_and_conc, -y) %>%
    mutate(columnName = gsub(fcStr, fcStrNorm, columnName) %>% 
             ## Preserve column order of the measurements:
             factor(levels = fcNormCols)) %>% 
    spread(columnName, yNormalized)
  
  ## ----------------------------------------------------------------------- ##
  ## Finalize
  ## ----------------------------------------------------------------------- ##
  
  ## Create output and preserve row names:
  out <- data %>% left_join(dataNormed, by = c(idVar, "temperature"))
  rownames(out) <- attr(data, "row.names")
  
  # Update attributes: Store column prefix of the new fold change columns
  importSettings$fcStrNorm = fcStrNorm
  attr(out, "importSettings") <- importSettings
  attr(out, "configTable")    <- attr(data, "configTable")
  
  message("Done.")
  return(out)
  
}
