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
  
  if (!missing(configTable)){
    warning("`configTable` is deprecated.", call. = TRUE)
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
  fcStr     <- importSettings$fcStr
  
  # Check whether intensitStr has class character.
  if (is.null(fcStr)){
    stop("attr(data, 'fcStr') must contain a field named 'fcStr'.")
  } else if (!is.character(fcStr)){
    stop("attr(data, 'importSettings')$fcStr must be of class character.")
  } else {
    message("Looking for fold change column prefix: '", fcStr, "'")
  }
  
  if (!any(grepl(fcStr, colnames(data)))){
    stop("The given prefix for fold change columns (specified by attr(data, 'importSettings')$fcStr) was not found in the column names of 'data'.")
  }
  
  fcStrNorm <- paste("norm", fcStr, sep="_")
  
  message("Performing median normalization...")
  
  norm.table <- do.call(rbind, lapply(unique(data$temperature), function(temp){
    
    # subset data to one temperature
    
    sub.table <- data[which(data$temperature==temp),]
    
    norm.list <- sapply(colnames(sub.table),function(coln){
      if (grepl(fcStr, coln)){
        col.median <- median(as.numeric(as.character(sub.table[[coln]])), 
                             na.rm=TRUE)
        norm.col <- as.numeric(as.character(sub.table[[coln]])) / col.median
        return(norm.col)
      }
    })
    
    nonNullFields <- norm.list[!sapply(norm.list, is.null)]
    
    # Create data frame from list
    # Caution: Set check.names = FALSE. Otherwise, special characters in column 
    # names are converted to '_' without warning -> can be problematic if
    # concentrations were imported from config file in scientific format. 
    # I.e. colname 'rel_fc_7.0000000000000007E-2' becomes 'rel_fc_7.0000000000000007E_2'
    norm.df <- data.frame(nonNullFields, check.names = FALSE) 
    
    return(norm.df)
  }))
  
  colnames(norm.table) <- gsub(fcStr, fcStrNorm, colnames(norm.table))
  
  out <- cbind(data, norm.table)
  
  # Update attributes: Store column prefix of the new fold change columns
  importSettings$fcStrNorm = fcStrNorm
  
  attr(out, "importSettings") <- importSettings
  
  attr(out, "configTable")    <- configTable
  
  message("Done.")
  return(out)
  
}
