#' @title Deprecated functions in package \sQuote{TPP}
#' @description These functions are provided for compatibility with older versions
#' of \sQuote{TPP} only, and will be defunct at the next release.
#'  
#' @details The following function is deprecated and will be made defunct; use
#'  the replacement indicated below:
#'  \itemize{
#'  \item{tpp2dPlotCCRGoodCurves: \code{\link{tpp2dCreateDRplots}}}
#'  \item{tpp2dPlotCCRSingleCurves: \code{\link{tpp2dCreateDRplots}}}
#'  \item{tpp2dPlotCCRAllCurves: \code{\link{tpp2dCreateDRplots}}}
#'  }
#' @export
#' @param configTable DEPRECATED
#' @param data DEPRECATED
#' @param idVar DEPRECATED 
#' @param fcStr DEPRECATED
#' @param verbose DEPRECATED
#' @param data.list DEPRECATED
#' @param intensityStr DEPRECATED
#' @return No value returned
#' 
#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dPlotCCRGoodCurves <- function(configTable=NULL, data=NULL, idVar="gene_name",
                                   fcStr="rel_fc_",  verbose=FALSE){
  
  .Deprecated("tpp2dCreateDRplots")
  
  plotList <- helperFctPlotGood(configTable = configTable, 
                                dataTable = data, 
                                idVar = idVar, 
                                fcStr = fcStr,
                                verbose = verbose,
                                paletteName = "Spectral")
  return(plotList)
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dPlotCCRAllCurves <- function(configTable=NULL, data=NULL, idVar="gene_name",
                                  fcStr="rel_fc_", verbose=FALSE){
  .Deprecated("tpp2dCreateDRplots")
  
  plotList <- helperFctPlotAll(configTable = configTable, 
                               dataTable = data, 
                               idVar = idVar, 
                               fcStr = fcStr,
                               verbose = verbose,
                               paletteName = "Spectral")
  
  return(plotList)
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dPlotCCRSingleCurves <- function(configTable=NULL, data=NULL, idVar="gene_name",
                                     fcStr="rel_fc_", verbose=FALSE){
  .Deprecated("tpp2dCreateDRplots")
  
  plotList <- helperFctPlotSingle(configTable = configTable, 
                                  dataTable = data, 
                                  idVar = idVar, 
                                  fcStr = fcStr,
                                  verbose = verbose)
  return(plotList)
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dEvalConfigTable <- function(configTable){
  # @title Evaluation of 2D-TPP Configuration File
  # @description Evaluates whether the configuration file is handed over as data frame or as file path 
  #   and loads the file path if necessary
  #  
  # @return A configtable that works with the 2D-TPP workflow
  # 
  # @examples 
  # # Import from data frame:
  # data(panobinostat_2DTPP_smallExample)
  # configTable <- tpp2dEvalConfigTable(panobinostat_2DTPP_config)
  # 
  # # Import from text file:
  # configPath <- system.file("test_data/panobinostat_ex_confg.txt", package = "TPP")
  # configTable <- tpp2dEvalConfigTable(configPath)
  # 
  # @param configTable data frame or character object with the path to a file, 
  #   that specifies important details of the 2D-TPP experiment. See Section 
  #   \code{details} for instructions how to create this object
  
  .Deprecated()
  
  # @title Evaluation of 2D-TPP Configuration File
  # @description Evaluates whether the configuration file is handed over as data frame or as file path 
  #   and loads the file path if necessary
  #  
  # @return A configtable that works with the 2D-TPP workflow
  # 
  # @param configTable data frame or character object with the path to a file, 
  #   that specifies important details of the 2D-TPP experiment. See Section 
  #   \code{details} for instructions how to create this object
  
  checkFunctionArgs(match.call(), c("configTable"))
  cfg <- importCheckConfigTable(infoTable = configTable, type = "2D")
  
  return(cfg)
  
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dRemoveZeroSias <- function(configTable, data.list, intensityStr="signal_sum_"){
  
  .Deprecated()
  
  data.list <- importFct_rmZeroSias(configTable, data.list, intensityStr)
  
  return(data.list)
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dReplaceColNames <- function(configTable, data.list, intensityStr, fcStr){
  
  .Deprecated()
  
  data.list <- importFct_createCCRInputFrom2DData(configTable, 
                                                  data.list, 
                                                  intensityStr, 
                                                  fcStr)
  return(data.list)
}

#' @name TPP-deprecated
#' @aliases TPP-deprecated
tpp2dCreateCCRConfigFile <- function(configTable){
  
  .Deprecated()
  
  out <- convert_2D_cfgTable_to_CCR_cfgTable(configTable = configTable)
  
  return(out) 
}
