#' @title Thermal proteome profiling (TPP)
#'   
#' @description \cite{TPP} is a toolbox for analyzing thermal proteome profiling
#'   (TPP) experiments.
#'   
#' @docType package
#' @name TPP
#' @param libname a character string giving the library directory where the 
#'   package defining the namespace was found. Passed to .onLoad function.
#' @param pkgname a character string giving the name of the package. Passed to 
#'   .onLoad function.
#'   
#' @references Savitski, M. M., Reinhard, F. B., Franken, H., Werner, T., 
#'   Savitski, M. F., Eberhard, D., ... & Drewes, G. (2014). Tracking cancer 
#'   drugs in living cells by thermal profiling of the proteome. Science, 
#'   346(6205), 1255784.
#'   
#'   Franken, H, Mathieson, T, Childs, D. Sweetman, G. Werner, T. Huber, W. & Savitski, M. M. (2015),
#'   Thermal proteome profiling for unbiased identification of drug targets and detection of downstream effectors.
#'   Nature protocols 10(10), 1567-1593.
#'   
#' @return No return value defined for this document.
#'   
#' @details In order to start a TPP-TR analysis, use function 
#'   \code{\link{analyzeTPPTR}}. For a TPP-CCR analysis, use function 
#'   \code{\link{analyzeTPPCCR}}. See the vignette for detailed instructions.

.onLoad <- function(libname, pkgname) {
  if (.Platform$OS.type == "windows") {
    if (Sys.which("zip")==""){
      msgText <- "\n ==> PLEASE READ BEFORE PACKAGE USE: 
     Could not locate a zip command in your path. This command is required for 
     Excel output on Windows (see vignette). 
     You can still use the 'TPP' package and access its results via the 
     dataframes produced by 'analyzeTPPTR' and 'analzyeTPPCCR'. 
     For export to Excel, please add a zip command to the path."
     packageStartupMessage(msgText,"\n")      
    }
  }
}
