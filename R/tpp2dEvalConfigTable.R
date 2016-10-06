#' @title Evaluation of 2D-TPP Configuration File
#' @description Evaluates whether the configuration file is handed over as data frame or as file path 
#'   and loads the file path if necessary
#'  
#' @return A configtable that works with the 2D-TPP workflow
#' 
#' @examples 
#' data("panobinostat_2DTPP_smallExample")
#' configTable <- tpp2dEvalConfigTable(panobinostat_2DTPP_config)
#' 
#' configTable <- tpp2dEvalConfigTable(system.file("example_data/2D_example_data/hdac2D_config.txt", 
#'                                                 package = "TPP"))
#' 
#' @param configTable data frame or character object with the path to a file, 
#'   that specifies important details of the 2D-TPP experiment. See Section 
#'   \code{details} for instructions how to create this object
#'   
#' @export
tpp2dEvalConfigTable <- function(configTable){
  cfg <- importCheckConfigTable(infoTable=configTable, type="2D")
  return(cfg)
}