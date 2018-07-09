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
#'   The default settings are adjusted to analyze data of the python package 
#'   \code{isobarQuant}. You can also customize them for your own dataset.
#'   
#'   The \code{configTable} argument is a dataframe, or the path to a 
#'   spreadsheet (tab-delimited text-file without quoted strings, or xlsx format). 
#'   Information about each experiment is stored row-wise. 
#'   It contains the following columns: 
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
#' be accessed by \code{Biobase::exprs(S)}. Protein expNames can be accessed by 
#' \code{featureNames(S)}. Isobaric labels and the corresponding concentrations are 
#'   returned by \code{S$label} and \code{S$concentration}.
#' 
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config, 
#' data = hdacCCR_data)
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

tppccrImport <- function(configTable, data=NULL, idVar="gene_name", 
                         fcStr="rel_fc_", naStrs=c("NA", "n/d", "NaN", "<NA>"), 
                         qualColName="qupm", nonZeroCols="qssm"){
  
  dataList <- importTR_main(configTable=configTable, data=data, idVar=idVar, 
                            fcStr=fcStr, naStrs=naStrs, qualColName=qualColName,
                            type="CCR")
  
  ## Remove proteins where nonZeroCols == 0:
  if (!is.null(nonZeroCols)){
    dataListFiltered <- list()
    for (expName in names(dataList)){
      message("Filtering CCR dataset: ", expName)
      
      dTmp <- dataList[[expName]]
      fDat <- pData(featureData(dTmp))
      
      jCol <- match(nonZeroCols, colnames(fDat))
      colStr <- paste("'", paste(nonZeroCols, collapse="', '"), "'", sep="")
      if (!any(is.na(jCol))){
        rm <- apply(as.matrix(fDat[, jCol]) == 0, 1, any)
        if(length(rm) > 0) {
          dTmpNew <- dTmp[!rm,]    
        }
        message("Removed proteins with zero values in column(s) ", colStr, ":")
        message("\t", nrow(dTmpNew), " out of ", nrow(dTmp), " proteins remaining.")
      } else {
        dTmpNew <- dTmp
        warning("At least one of the column(s) ", colStr, 
                ", which were specified in the argument 'nonZeroCols' could not be found in the data. No filtering performed.")
      }
      
      dataListFiltered[[expName]] <- dTmpNew
    }
  } else {
    dataListFiltered <- dataList
  }

  
  ## Convert given concentrations to log scale for curve fitting and plotting:
  for (expName in names(dataListFiltered)){
    datTmp <- dataListFiltered[[expName]]
    concTmp <- log10(as.numeric(as.character(datTmp$concentration)) * 10^-6)
    concTmp[is.infinite(concTmp)] <- -15  
    dataListFiltered[[expName]]$concentration <- concTmp
  }
  
  
  return(dataListFiltered)
}