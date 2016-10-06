#' @title Remove rows with zero sumionarea values 
#' 
#' @description Removes zero sumionare values in a specified data.list so that no errors are
#'   generated in the following fold change computation step. A corresponding data.list  
#'   with NAs instead of zeros is returned.
#'   
#' @return A list of data frames with NAs instead of zeros.
#'
#' @param configTable data frame that specifies important details of the 2D-TPP experiment.
#' @param data.list list of data frames of corresponding experiment data  
#' @param intensityStr character string indicating which columns contain the sumionarea 
#'   values. Those column names containing the suffix \code{intensityStr} 
#'   will be regarded as containing sumionare values.
#'   
#' @export 
tpp2dRemoveZeroSias <- function(configTable, data.list, intensityStr="sumionarea_protein_"){
  lapply(names(data.list), function(l.name){
    # get sumionare cloumns
    intensity.cols <- colnames(data.list[[l.name]])[grep(intensityStr, colnames(data.list[[l.name]]))]
    new.intensity.df <- data.frame(sapply(intensity.cols, function(scol){
      num.vals <- as.numeric(sapply(data.list[[l.name]][scol], function(val) return(as.character(val))))
      zero.ids <- which(num.vals==0)
      if(!length(zero.ids)==0){
        num.vals[zero.ids] <- NA
        return(num.vals)
      }else{
       return(num.vals) 
      }
    }))
    data.list[[l.name]][intensity.cols] <- new.intensity.df # new
  })
  return(data.list)
}