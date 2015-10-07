#' @title Normalize protein fold changes
#'   
#' @description Normalizes fold changes determined by TPP-TR experiments over 
#'   different experimental groups.
#'   
#' @return A list of ExpressionSets storing the normalized data for each
#'   experiment. Each ExpressionSet contains the measured fold changes, as well
#'   as row and column metadata. In each ExpressionSet \code{S}, the fold
#'   changes can be accessed by \code{exprs(S)}. Protein expNames can be
#'   accessed by \code{featureNames(S)}. Isobaric labels and the corresponding temperatures are 
#'   returned by \code{S$label} and \code{S$temperature}
#'   
#' @details Performs normalization of all fold changes in a given list of 
#'   ExpressionSets. The normalization procedure is described in detail in 
#'   Savitski et al. (2014). Whether normalization needs to be performed and 
#'   what method is best suited depends on the experiment. Here we provide a 
#'   reasonable solution for the data at hand.
#'   
#'   We distinguish between filtering conditions on fold changes and on 
#'   additional annotation columns. Correspondingly, \code{normReqs} contains 
#'   two fields, \code{fcFilters} and \code{otherFilters}. Each entry contains a
#'   data frame with three columns specifying the column to be filtered, as well
#'   as upper and lower bounds. An example is given by 
#'   \code{\link{tpptrDefaultNormReqs}}.
#'   
#' @param data List of \code{ExpressionSet}s with protein fold changes to be 
#'   normalized.
#' @param normReqs List of filtering criteria for construction of the 
#'   normalization set.
#' @param qcPlotTheme ggplot theme for the created plots
#' @param qcPlotPath location where plots of the curves fitted to the 
#'   normalization set medians should be stored.
#' @param startPars start values for the melting curve parameters. Will be 
#'   passed to function \code{\link{nls}} for curve fitting.
#' @param maxAttempts maximal number of curve attempts to fit melting curve to 
#'   fold change medians when computing normalization factors.
#' @param fixedReference name of a fixed reference experiment for normaliztion. 
#'   If NULL (default), the experiment with the best R2 when fitting a melting 
#'   curve through the median fold changes is chosen as the reference.
#' @seealso \code{\link{tpptrImport}}
#'   
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=tpptrDefaultNormReqs())
#' names(tpptrNorm)
#' 
#' @references Savitski, M. M. and Reinhard, F. BM. and Franken, H. and Werner, 
#'   T. and Savitski, M. F. and Eberhard, D. and Molina, D. M. and Jafari, R. 
#'   and Dovega, R. B. and Klaeger, S. and others (2014) Tracking cancer drugs 
#'   in living cells by thermal profiling of the proteome. Science 346(6205), p.
#'   1255784.
#'   
#' @export
tpptrNormalize <- function(data, normReqs=tpptrDefaultNormReqs(), 
                      qcPlotTheme=tppDefaultTheme(), qcPlotPath=NULL, 
                      startPars=c("Pl"=0, "a"=550, "b"=10), maxAttempts=1, fixedReference=NULL){
  ## 1. Filter data and detect treatment group with most remaining proteins:
  infoNormP <- filterTables(data=data, normReqs=normReqs)
  
  ## 2. Reduce all data sets to only those proteins contained in normP:
  normP     <- infoNormP[["protein_IDs"]]
  listNormP <- sapply(data, function(x){exprSubset(x, subset=normP)}, simplify=FALSE)
  message("-----------------------------------")
  ## 3. Fit sigmoids to medians of each treatment group and determine best fit:
  listNormFit <- computeNormFactors(data=listNormP, startPars=startPars, maxAttempts=maxAttempts, fixedReference=fixedReference)
  
  ## 4. Create QC plot of melting curve fit
  meltCurveModels <- listNormFit[["models"]]
  fcMedians       <- listNormFit[["medians"]]
  tempVals        <- listNormFit[["tempVals"]]
  r2Vec           <- listNormFit[["rSquared"]]
  nNormP          <- length(normP)
  
  qcPlotsMedianFits <- plotNormCurves(modelList=meltCurveModels, xMat=tempVals, 
                                      fcMat=fcMedians, r2Vec=r2Vec, 
                                      nNormP=nNormP, plotTheme=qcPlotTheme)
  if (!is.null(qcPlotPath)){
    pdf(file=file.path(qcPlotPath, "QCplots_median_fits.pdf"), width=8, height=9)
    print(qcPlotsMedianFits)
    dev.off()
  }
  
  message("-----------------------------------")
  ## 5. Normalize all proteins by best curve:
  message("Normalizing all proteins in all experiments.")
  dfCoeffs = listNormFit[["corrFactors"]]
  normData <- sapply(names(data), simplify=FALSE,
                     function(n){applyCoeffs(data=data[[n]], coeffs=dfCoeffs[,n])})
  
  ## 6. Return results:
  message("Normalization successfully completed!\n")
  return(list(normData=normData, qcPlotObj=qcPlotsMedianFits))
}