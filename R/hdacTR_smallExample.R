#' @title Example subset of a Panobinostat TPP-TR dataset and the
#'   corresponding configuration table to start the analysis.
#' @name hdacTR_smallExample
#' @docType data
NULL
#' 
#' @rdname hdacTR_data
#' @name hdacTR_data
#' @title TPP-TR example dataset.
#' @description Example subset of a dataset obtained by TPP-TR experiments to 
#'   investigate possible targets for HDAC inhibitor Panobinostat.
#' @details \code{hdacTR_data} is a list of data frames that contain
#'   measurements for HDACs as well as a random selection of 500 further
#'   proteins.
#'   
#'   You can use this dataset to explore the \code{\link{TPP}} package 
#'   functionalities without invoking the whole time consuming analysis on the 
#'   whole dataset.
#'   
#'   The originial dataset is located in the folder 
#'   \code{'example_data/TR_example_data'} in the package's installation
#'   directory. You can find it on your system by the \code{R} command
#'   \code{system.file('example_data', package = 'TPP')}.
NULL
#'
#'@rdname hdacTR_config
#'@name hdacTR_config
#'@title The configuration table to analyze \link{hdacTR_data}.
#'@description The configuration table to analyze \link{hdacTR_data}.
#'@details \code{hdacTR_config} is a data frame that specifies the experiment 
#'  name, isobaric labels, and the administered temperatures at each 
#'  label.
NULL 