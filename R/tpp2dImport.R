#' @title Import 2D-TPP data
#' 
#' @description Imports data from 2D-TPP experiments by parsing a configTable and reading in 
#'   corresponding data file or data frames containing raw data (sumionarea values) and creating a 
#'   big data frame comprising all samples with respective fold changes
#' 
#' @return A dataframe comprising all experimental data
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
#' @param nonZeroCols character string indicating a column that will be used for
#'   filtering out zero values.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param addCol additional column names that specify columns in the input data that are 
#'   to be attached to the data frame throughout the analysis 
#' 
#' @export
tpp2dImport <- function(configTable = NULL, 
                        data = NULL, 
                        idVar = "gene_name",
                        addCol = NULL,
                        intensityStr = "signal_sum_",   
                        qualColName = "qupm", 
                        nonZeroCols = "qssm",
                        fcStr = NULL){
  
  # check configTable and read in table if necessary
  configWide <- importCheckConfigTable(infoTable = configTable, type = "2D")
  #configLong <- gatherConfigTable(config = configWide)
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  experiment = unique_ID <- NULL
  
  message("Importing data...")
  # import data as list, if is.null(fcStr) function will omit fold changes and only read in
  # sumionareas
  dataList <- import2DTR_main(configTable = configWide, 
                              data = data, 
                              idVar = idVar, 
                              fcStr = fcStr,
                              addCol = addCol,
                              naStrs = c("NA", "n/d", "NaN"),
                              intensityStr = intensityStr,
                              qualColName = qualColName,
                              nonZeroCols = nonZeroCols)
  
  # create one wide data table
  dataWide <- importFct_createCCRInputFrom2DData(configTable = configWide, 
                                                  data.list = dataList, 
                                                  intensityStr = intensityStr, 
                                                  fcStr = fcStr) %>%
    mutate(experiment = factor(experiment), unique_ID = factor(unique_ID))
  
  dataWide <- dataWide %>% arrange_(.dots = c(idVar, "temperature"))
  
  ## Add annotation for use in later functions:
  attr(dataWide, "configTable") <- configWide
  
  importSettings <- list(proteinIdCol = idVar, 
                         uniqueIdCol = "unique_ID",
                         addCol = addCol,
                         intensityStr = intensityStr, 
                         qualColName = qualColName, 
                         nonZeroCols = nonZeroCols,
                         fcStr = fcStr)
  attr(dataWide, "importSettings") <- importSettings
  
  message("Done.")
  return(dataWide)
}
