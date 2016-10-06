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
#'                                               dataTable = headData2d, 
#'                                               intensityStr = "sumionarea_protein_")  
#'                                  
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param dataTable dataframe that contain the data for the 2D-TPP experiment
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the suffix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param fcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{fcStr} will be pasted in front of
#'   the names of the experiments.
#' 
#'   
#' @export
tpp2dComputeFoldChanges <- function(configTable=NULL, dataTable=NULL, intensityStr=NULL, 
                                 fcStr="rel_fc_protein_"){
  if (!is.null(intensityStr) && is.character(intensityStr)){
    message("Computing fold changes...")
    # determine reference colnames for experiments
    intensity.col.ids <- grep(intensityStr, colnames(dataTable))
    fc.cols <- sapply(intensity.col.ids, function(sc){
      return(paste(fcStr, sub(intensityStr, "", colnames(dataTable)[sc]), sep=""))
    })
    ref.col.id <- grep(paste(intensityStr, "0$", sep=""), colnames(dataTable))
    if (length(ref.col.id)==0){
      ref.col.id <- grep(paste(intensityStr, "0.0$", sep=""), colnames(dataTable))
      if (length(ref.col.id)==0){
        ref.col.id <- grep(paste(intensityStr, "0.00$", sep=""), colnames(dataTable))
        if (length(ref.col.id)==0){
          ref.col.id <- grep(paste(intensityStr, "0.000$", sep=""), colnames(dataTable))
        }
      }
    }
    ref.vals <- as.numeric(as.character(dataTable[,ref.col.id])) 
    # calc fold change of respective values by dividing by the reference value
    fc.m <- data.frame(matrix(as.numeric(as.character(unlist(dataTable[,intensity.col.ids]))),
                              nrow=nrow(dataTable[,intensity.col.ids])) / ref.vals)
    names(fc.m) <- fc.cols
    dataTable[fc.cols] <- fc.m
    message("Done.")
    return(dataTable)
  }else{
    stop("You must specifiy a valid intensityStr!")
  }
} 
