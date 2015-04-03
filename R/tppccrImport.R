#' @title Import TPP-CCR dataset for analysis by the 
#'   \code{\link{TPP}} package.
#'   
#' @description \code{tppccrImport} imports a table of protein fold changes and 
#' stores them in an ExpressionSet for use in the \code{\link{TPP}} package.
#'   
#' @details The imported dataset has to contain measurements obtained by a
#'   TPP-CCR experiment. Fold changes need to be pre-computed using the lowest 
#'   concentration as reference.
#'   
#'   The dataset can be specified by filename in the \code{configTable} 
#'   argument, or given directly in the \code{data} argument
#'   
#'   The default settings are adjusted to analyse data of the pyhton package 
#'   \code{isobarQuant}. You can also customise them for your own dataset.
#'   
#'   The \code{configTable} argument is a dataframe, or the path to a 
#'   spreadsheet (tab-delimited text-file or xlsx format). Information about 
#'   each experiment is stored row-wise. It contains the following columns: 
#'   \itemize{
#'  \item{\code{Path}: }{location of the datafile. Alternatively, data can be directly handed
#'   over by the \code{data} argument.}
#'  \item{\code{Experiment}: }{unique experiment name.}
#'  \item{Label columns: } each isobaric label names a column that contains the 
#'   concentration administered for the label in the individual experiments.
#' }
#'   
#'   During data import, proteins with NAs in the data column specified by \code{idVar} receive 
#'   unique generic IDs so that they can be processed by the package.
#'
#' @return ExpressionSet object storing the measured fold changes, as well as
#' row and column metadata. In each ExpressionSet \code{S}, the fold changes can 
#' be accessed by \code{exprs(S)}. Protein expNames can be accessed by 
#' \code{featureNames(S)}. TMT labels and the corresponding temperatures are 
#' returned by \code{S$labels} and \code{S$temperatures}.   
#' 
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config_repl1, data = hdacCCR_data_repl1)
#' 
#' @param configTable either a dataframe or the path to a spreadsheet. In both cases
#' it specifies necessary information of the TPP-CCR experiment.
#' @param data dataframe containing fold change measurements and 
#' additional annotation columns to be imported. Can be used instead of 
#' specifying the file path in \code{configTable}.
#' @param idVar character string indicating which data column provides the unique 
#' identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix 
#'   \code{fcStr} will be regarded as containing fold change values.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param nonZeroCols character string indicating a column that will be used for
#' filtering out zero values.
#' @export
#' @seealso \code{\link{tpptrImport}}, \code{\link{tppccrCurveFit}}

tppccrImport <- function(configTable, data=NULL, idVar="gene_name", fcStr="rel_fc_", 
                       naStrs=c("NA", "n/d", "NaN", "<NA>"), qualColName="qupm", nonZeroCols="qssm"){
  message("Importing data...\n")
  configTable <- importCheckConfigTable(infoTable=configTable)
  
  expName        <- configTable$Experiment
  file           <- configTable$Path
  labels         <- setdiff(colnames(configTable), c("Experiment", "Path"))
  concentrations <- configTable[, labels]
  
  ## Determine dataset name:
  expNameCheck   <- importCheckExperimentNames(expNames=expName, expectedLength=1)
  expName        <- expNameCheck$expNames
  genericExpNames <- expNameCheck$genericExpNames
  
  ## Check whether dataframe or filename is specified for data import:
  argList <- importCheckDataFormat(dataframes=data, files=file, 
                                   expNames=expName, genericExpNames=genericExpNames)
  data <- argList[["dataframes"]][[expName]]
  file      <- argList[["files"]][[expName]]
  
  ## Determine matrix of concentrations to each isotope label
  concMatrix <- importCheckTemperatures(temp=concentrations, nData=1)
  concVec    <- as.numeric(concMatrix)
  
  ## Check isotope label argument
  labels <- importCheckLabels(labelValues=labels, temperatureMat=concMatrix)
  
  ## Import data and convert into ExpressionSet format:
  eSet <- importSingleExp(name=expName, filename=file, dataframe=data,
                          labels=labels, labelValues=concVec, type="CCR",
                          idVar=idVar, fcStr=fcStr, qualColName=qualColName, naStrs=naStrs)
  
  ## ---------------------------------------------------------------------------
  ## Preprocess imported dataset
  ## ---------------------------------------------------------------------------
  ## Filter data where qusm = 0:
  filterCol <- match(nonZeroCols, varLabels(featureData(eSet)))
  if (!is.na(filterCol)){
    message("Removing proteins with zero values in column ", paste(varLabels(featureData(eSet))[filterCol], collapse=', '), ":")
    rm <- pData(featureData(eSet))[, filterCol] == 0
    if(length(rm) > 0) {
      eSetFiltered <- eSet[!rm,]    
    }
    message(nrow(eSetFiltered), " out of ", nrow(eSet), " proteins remaining after filtering.")  
    message("\n")    
  } else {
    eSetFiltered <- eSet
    warning("Column '", nonZeroCols, "', which was specified in the argument 'nonZeroCols' could not be found in the data. No filtering performed.")
  }
  return(eSetFiltered)
}