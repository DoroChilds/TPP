#' @title Example subsets of a Panobinostat TPP-CCR dataset (replicates 1 and 2) and the
#'   corresponding configuration table to start the analysis.
#' @name hdacCCR_smallExample
#' @docType data
NULL
#' 
#' @rdname hdacCCR_data
#' @name hdacCCR_data
#' @title TPP-CCR example dataset (replicates 1 and 2)
#' @description Example subset of a Panobinostat TPP-CCR dataset (replicates 1 
#' and 2)
#' @details A list with two subsets of a dataset obtained by TPP-CCR experiments to 
#'   investigate drug effects for HDAC inhibitor Panobinostat. It contains 7 
#'   HDACs as well as a random selection of 493 further proteins.
#'   
#'   You can use this dataset to explore the \code{\link{TPP}} package 
#'   functionalities without invoking the whole time consuming analysis on the 
#'   big dataset.
#'   
#'   The originial dataset is located in the folder 
#'   \code{'example_data/CCR_example_data'} in the package's installation 
#'   directory. You can find it on your system by the \code{R} command 
#'   \code{system.file('example_data', package = 'TPP')}.
NULL
#' @rdname hdacCCR_config
#' @name hdacCCR_config
#' @title The configuration table to analyze \link{hdacCCR_data}.
#' @description The configuration table to analyze \link{hdacCCR_data}.
#' @details \code{hdacCCR_config} is a data frame that specifies the experiment
#'  names, isobaric labels, and the administered drug concentrations at each
#'  label.
NULL