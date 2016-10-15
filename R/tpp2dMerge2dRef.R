#' @title Merge 2D-TPP result data with TPP-TR reference data
#'   
#' @description Merges 2D-TPP result data with TPP-TR reference data to generate
#' a big table including both results
#'  
#' @return A data frame with results merged from 2D-TPP and TPP-TR reference 
#' 
#' @param data dataframe containing the 2D-TPP results
#' @param trRef character string of a valid system path to a TPP-TR reference 
#' RData object
#' @param idVar character string matching the column containing the unique protein 
#' identifiers
#' 
tpp2dMerge2dRef <- function(data=NULL, trRef=NULL, idVar="representative"){
  if (!is.null(data) && !is.null(trRef)){
    # load trRef
    load(trRef)
    sumResTable <- subset(tppRefData$sumResTable$summary, select=-c(clustername))
    
    # merge 2D with trRef
    data <- merge(data, sumResTable, by.x=idVar, by.y="Protein_ID", 
                            all.x=TRUE, all.y=FALSE)
    
    # return resulting table
    return(data)
    
  }else {
    stop("Please specify a valid argument for 'data' and 'trRef'!")
  }
}