#' @title Import 2D-TPP data
#' 
#' @description Imports data from 2D-TPP experiments by parsing a configTable and reading in 
#'   corresponding data file or data frames containing raw data (sumionarea values) and creating a 
#'   big data frame comprising all samples with respective fold changes
#' 
#' @return A dataframe comprising all experimental data
#' 
#' @examples 
#'   data("panobinostat_2DTPP_smallExample")
#'   tpp2dResults <- tpp2dImportData(configTable = panobinostat_2DTPP_config, 
#'                              data = panobinostat_2DTPP_data, fcStr = NULL)
#' 
#' @param configTable dataframe, or character object with the path to a file, 
#'   that specifies important details of the 2D-TPP experiment. See Section 
#'   \code{details} for instructions how to create this object.
#' @param data single dataframe, containing raw measurements and if already available fold
#'   changes and additional annotation columns to be imported. Can be used instead of 
#'   specifying the file path in the \code{configTable} argument.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the suffix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param addCol additional column names that specify columns in the input data that are 
#'   to be attached to the data frame throughout the analysis 
#' 
#' @export
tpp2dImportData <- function(configTable=NULL, data=NULL, idVar="representative", 
                            addCol=c("clustername", "msexperiment_id"),
                            intensityStr="sumionarea_protein_", qualColName=c("qupm","qusm"),
                            fcStr="rel_fc_protein_"){
  if (is.null(configTable) || !is.data.frame(configTable)){
    stop("Please specify a valid configTable of type data.frame!")
  }
  message("Importing data...")
  # import data as list, if is.null(fcStr) function will omit fold changes and only read in
  # sumionareas
  Data2d <- tpp2dCreateDataFrameList(configTable=configTable, data=data, 
                                     idVar=idVar, fcStr=fcStr, addCol=addCol,
                                     intensityStr=intensityStr, qualColName=qualColName)
  message("Done.")
  # remove 0 sumionarea values
  d.list <- tpp2dRemoveZeroSias(configTable=configTable, data.list=Data2d, 
                    intensityStr=intensityStr)
  
  # create one big data table
  d.table <- tpp2dReplaceColNames(configTable=configTable, data.list=d.list, 
                                  intensityStr=intensityStr, fcStr=fcStr)
  message("Done.")
  return(d.table)
}