#' @title Import TPP-TR datasets for analysis by the \code{\link{TPP}} package.
#'   
#' @description \code{tpptrImport} imports several tables of protein fold 
#'   changes and stores them in a list of ExpressionSets for use in the 
#'   \code{\link{TPP}} package.
#'   
#' @details The imported datasets have to contain measurements obtained by 
#'   TPP-TR experiments. Fold changes need to be pre-computed using the lowest 
#'   temperature as reference.
#'   
#'   An arbitrary number of datasets can be specified by filename in the 
#'   \code{Path}-column of the \code{configTable} argument, or given directly as
#'   a list of dataframes in the \code{data} argument. They can differ, for 
#'   example, by biological replicate or by experimental condition (for example,
#'   treatment versus vehicle). Their names are defined uniquely by the 
#'   \code{Experiment} column in \code{configTable}. Experimental conditions can
#'   be specified by an optional column in \code{configTable}.
#'   
#'   The default settings are adjusted to analyze data of the python package 
#'   \code{isobarQuant}. You can also customize them for your own dataset.
#'   
#'   The \code{configTable} argument is a dataframe, or the path to a 
#'   spreadsheet (tab-delimited text-file without quoted strings, or xlsx format). 
#'   Information about each experiment is stored row-wise. 
#'   It contains the following columns: 
#'   \itemize{ \item{\code{Path}:}{location of each datafile. Alternatively, 
#'   data can be directly handed over by the \code{data} argument.} 
#'   \item{\code{Experiment}: }{unique experiment names.} 
#'   \item{\code{Condition}: }{experimental conditions of each dataset.} 
#'   \item{Label columns: } each isobaric label names a column that contains the
#'   temperatures administered for the label in the individual experiments. }
#'   
#'   Proteins with NAs in the data column specified by \code{idVar} receive 
#'   unique generic IDs so that they can be processed by the package.
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(hdacTR_config, hdacTR_data)
#'   
#' @param configTable either a dataframe or the path to a spreadsheet. In both 
#'   cases it specifies necessary information of the TPP-CCR experiment.
#' @param data single dataframe, or list of dataframes, containing fold change 
#'   measurements and additional annotation columns to be imported. Can be used 
#'   instead of specifying the file path in \code{configTable}.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param outputFormat output format. Either "eSetList" to obtain output in the 
#'  same way as previously (will be deprecated soon), or "tidy" to obtain a 
#   tidy table of fold changes (the recommended setting).

#' @return A list of ExpressionSets storing the imported data for experiment.
#'   Each ExpressionSet contains the measured fold changes, as well as row and
#'   column metadata. In each ExpressionSet \code{S}, the fold changes can be
#'   accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#'   \code{featureNames(S)}. Isobaric labels and the corresponding temperatures are 
#'   returned by \code{S$label} and \code{S$temperature}
#'   
#' @export
#' @seealso \code{\link{tppccrImport}}

tpptrImport <- function(configTable, data=NULL, idVar="gene_name", 
                        fcStr="rel_fc_", naStrs=c("NA", "n/d", "NaN"), 
                        qualColName="qupm", outputFormat = "eSetList"){
  if (outputFormat == "eSetList"){
    # warning("The outputFormat 'eSetList' is deprecated and will be removed soon.\n  Use outputFormat='tidy' instead.")
    out <- importTR_main(configTable=configTable, data=data, idVar=idVar, 
                         fcStr=fcStr, naStrs=naStrs, qualColName=qualColName,
                         type="TR")
    
  } else if (outputFormat == "tidy"){
    out <- importTR_tidy(configTable = configTable, data = data, idVar = idVar, 
                         fcStr = fcStr, naStrs = naStrs, qualColName = qualColName,
                         type = "TR")
  }
  return(out)
}
