#' @title Create Report of 2D-TPP analysis
#' 
#' @description Creates a markdown pdf file that summarizes the 2D-TPP analysis by reporting e.g. R
#'  version and package versions used
#' 
#' @return A pdf or html report which summarizes all parameters that were set
#' 
#' @param resultPath character string containing a system path to where the report should be 
#'  written
#' @param configFile character string containing a valid system path to a file which summarizes 
#'  the experimental details of the 2D-TPP experiment or respective data frame
#' @param configTable data frame summarizing the experimental details of the 2D-TPP experiment
#' @param resultTable output data frame from an 2D-TPP analysis
#' @param idVar unique protein identifier prefix
#' @param fcStr fold change identifier prefix
#' @param intensityStr intensity values prefix
#' @param addCol vector of strings indicating which additional data columns were imported
#' @param fcTolerance tolerance for the fcCutoff parameter
#' @param r2Cutoff Quality criterion on dose response curve fit.
#' @param fcCutoff Cutoff for highest compound concentration fold change
#' @param slopeBounds Bounds on the slope parameter for dose response curve 
#'   fitting
#' @param fTest boolean variable stating whether an fTest was peformed 
#' @param trRef character string containing a valid system path to a previously generated TPP-TR
#'  reference object
#' @param documentType character string indicating which document type the report should have
#'  default: "html_document", alternatives: "pdf_document"
#' @param normalize boolean flag indicating whether median normalization has been performed 
#' @param methods vector of characters which indicate which methods have been used 
#' @param fcStrUpdated character string matching the fold change columns after normalization has
#'  been performed 
#' 
#' @export
tpp2dCreateReport <- function(resultPath=NULL, configFile=NULL, documentType="html_document",
                              configTable=NULL, normalize=TRUE, methods=c(""),
                              resultTable=NULL, idVar=NULL, fcStr = "rel_fc_protein",
                              fcStrUpdated = "norm_rel_fc_protein_", 
                              intensityStr=NULL, addCol=NULL,
                              fcTolerance=NA, r2Cutoff=NA, fcCutoff=NA, slopeBounds=c(NA,NA),
                              fTest=FALSE, trRef="none"){
  
  if(!is.null(resultPath) && file.exists(resultPath)){
    message("Creating report...\n")
    inFile <- file.path(system.file(package="TPP"), "tpp2d_report.Rmd")
    if (documentType=="html_document"){
      outFile <- file.path(resultPath, paste("report_", basename(resultPath), ".html", sep=""))
    }else if(documentType=="pdf_document"){
      outFile <- file.path(resultPath, paste("report_", basename(resultPath), ".pdf", sep=""))
    }else{
      stop("Please specify a valid documentType!")
    }
    render(input=inFile , output_file=outFile, envir=environment(), intermediates_dir=tempdir(), 
           output_format=documentType)
  }else {
    stop("Please specify a valid resultPath!")
  }
}
