#' @title Compute 2D-TPP fold changes
#' @description Computes fold changes by calculating fold changes of the sumionarea 
#'  relative to the reference column.
#'
#' @return A data.frame with addtional columns with constitute fold changes calculated with 
#'  respect to the intenstiy values of the zero treatment column 
#' 
#' @examples
#' data("panobinostat_2DTPP_smallExample")
#' load(system.file("example_data/2D_example_data/referenceReadInData.RData", package="TPP"))
#' exampleFoldChanges <- tpp2dComputeFoldChanges(configTable = panobinostat_2DTPP_config, 
#'                                               data = headData2d, 
#'                                               intensityStr = "sumionarea_protein_")  
#'                                  
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param data dataframe that contain the data for the 2D-TPP experiment
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the suffix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param fcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{fcStr} will be pasted in front of
#'   the names of the experiments.
#' 
#'   
#' @export
tpp2dComputeFoldChanges <- function(configTable=NULL, data=NULL, intensityStr=NULL, 
                                 fcStr="rel_fc_protein_"){
  if (!is.null(intensityStr) && is.character(intensityStr)){
    message("Computing fold changes...")
    # determine reference colnames for experiments
    intensity.col.ids <- grep(intensityStr, colnames(data))
    fc.cols <- sapply(intensity.col.ids, function(sc){
      return(paste(fcStr, sub(intensityStr, "", colnames(data)[sc]), sep=""))
    })
    ref.col.id <- grep(paste(intensityStr, "0$", sep=""), colnames(data))
    if (length(ref.col.id)==0){
      ref.col.id <- grep(paste(intensityStr, "0.0$", sep=""), colnames(data))
      if (length(ref.col.id)==0){
        ref.col.id <- grep(paste(intensityStr, "0.00$", sep=""), colnames(data))
        if (length(ref.col.id)==0){
          ref.col.id <- grep(paste(intensityStr, "0.000$", sep=""), colnames(data))
        }
      }
    }
    ref.vals <- as.numeric(as.character(data[,ref.col.id])) 
    # calc fold change of respective values by dividing by the reference value
    fc.m <- data.frame(matrix(as.numeric(as.character(unlist(data[,intensity.col.ids]))),
                              nrow=nrow(data[,intensity.col.ids])) / ref.vals)
    names(fc.m) <- fc.cols
    data[fc.cols] <- fc.m
    message("Done.")
    return(data)
  }else{
    stop("You must specifiy a valid intensityStr!")
  }
} 
