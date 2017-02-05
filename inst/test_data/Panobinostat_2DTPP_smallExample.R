#' @title Example subsets of a Panobinostat 2D-TPP dataset and the
#'   corresponding configuration table to start the analysis.
#' @name Panobinostat_2DTPP_smallExample
#' @docType data
#' @description Example dataset obtained by 2D-TPP experiments for analysis by 
#' the TPP-package. It contains all necessary arguments to start the analysis
#' (config table and list of data frames).

NULL
#' 
#' @rdname panobinostat_2DTPP_data
#' @name panobinostat_2DTPP_data
#' @title 2D-TPP-CCR example dataset 
#' @description Example subset of a Panobinostat 2D-TPP dataset 
#' @details A list with two subsets of a dataset obtained by 2D-TPP experiments to 
#'   investigate drug effects for HDAC inhibitor Panobinostat. The experiment 
#'   was performed on living HepG2 cells (see Becher et al. (2016). Thermal 
#'   profiling reveals phenylalanine hydroxylase as an off-target of panobinostat. 
#'   Nature Chemical Biology, (September)) 
#'   It contains 7 HDACs as well as a random selection of 493 further proteins.
#'   
#'   You can use this dataset to explore the \code{\link{TPP}} package 
#'   functionalities without invoking the whole time consuming analysis on the 
#'   big dataset.
#'   
#'   The original dataset in plain text format is located in the folder 
#'   \code{'example_data/2D_example_data'} in the package's installation 
#'   directory. You can find it on your system by the \code{R} command 
#'   \code{system.file('example_data', package = 'TPP')}.
NULL
#' @rdname panobinostat_2DTPP_config
#' @name panobinostat_2DTPP_config
#' @title The configuration table to analyze \link{panobinostat_2DTPP_data}.
#' @description The configuration table to analyze \link{panobinostat_2DTPP_data}.
#' @details \code{panobinostat_2DTPP_config} is a data frame that specifies the experiment
#'  names, isobaric labels, and the administered drug concentrations at each
#'  label.
NULL