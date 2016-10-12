#' @title Plot quality control pEC50 plots
#'   
#' @description Plots quality control plots which indicate at which temperatures the pEC50 
#'  values of the treatment curves lie in comparison to those of the reference data
#'  
#' @return A folder with plots for each identified protein that compare melting points in
#'  the reference data set with the 2D-TPP data set
#'  
#' @param resultTable data.frame containing the results of a CCR analysis of 2D-TPP data
#' @param resultPath character string containing a valid system path to which the the qc
#'  plots will be written
#' @param trRef character string with a link to a TPP-TR reference object RData file
#' @param idVar character string indicating how the column that contains the unique protein 
#'  identifiers is called
#' 
tpp2dPlotQCpEC50 <- function(resultTable=NULL, resultPath=NULL, trRef=NULL, 
                             idVar="representative"){
  message("Creating melting point vs. pEC50 QC plots...")  
  
  
  # load TR reference object
  tppTRObj <- tpp2dTRReferenceObject(tppRefDataPath=trRef)
  # create resultTable subset of only rows have valid pEC50
  subCCR <- resultTable[which(!is.na(resultTable$pEC50)),]
  
  # loop over all protein IDs and plot melting point against pEC50
  plots <- lapply(as.character(unique(subCCR[[idVar]])), function(pID){
    #prot <- as.character(unique(subCCR[which(subCCR[[idVar]]==pID),][[protName]]))
    temp.df <- data.frame(
      temp=as.numeric(sub("C", "", subCCR[which(subCCR[[idVar]]==pID),][["temperature"]])),
      pEC50=subCCR[which(subCCR[[idVar]]==pID),][["pEC50"]])
    protCCRData <- list(protID=pID,
                      efficacyData=temp.df)
    qcPlot <- try(tppTRObj$createMeltPpEC50plot(protCCRData))
    if (class(qcPlot)[1]!="try-error"){
      if (!is.null(resultPath)){
        dirPath <- file.path(resultPath, "qc_pEC50plots")
        if (!file.exists(dirPath)){
          dir.create(dirPath)
        }
        filePath <- file.path(dirPath, paste(gsub("\\.", "_", pID), "qc_pEC50.pdf", sep="_"))
        ggsave(plot=qcPlot, filename=filePath, width=9, height=9)
      }
      return(qcPlot)
    }
  })
  names(plots) <- as.character(unique(subCCR[[idVar]]))
  return(plots)
  message("Done.")
}