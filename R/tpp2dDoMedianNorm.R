#' @title Median normalization of protein fold changes of 2D-TPP data
#'   
#' @description Normalizes fold changes retrieved from 2D-TPP experiment by dividing by the median fold
#'   change 
#'   
#' @examples 
#'   load(system.file("example_data/2D_example_data/referenceFcData.RData", package="TPP"))
#'   data("panobinostat_2DTPP_smallExample")
#'   NormData2d <- tpp2dDoMedianNorm(configTable = panobinostat_2DTPP_config, 
#'                                     dataTable=headData2dFc,
#'                                     fcStr="rel_fc_protein_")
#'   
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param dataTable data frame that contains the data for the 2D-TPP experiment
#' @param fcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{fcStr} will be pasted in front of
#'   the names of the experiments.
#'   
#' @return A dataframe identical to the input dataframe except that the columns containing the
#'   fold change values have been normalized by their median.
#'   
#' @export 
tpp2dDoMedianNorm <- function(configTable, dataTable, fcStr="rel_fc_protein_"){
  if (!any(grepl(fcStr, colnames(dataTable)))){
    stop("Please specify a valid fcStr suffix matching the fold change columns!")
  } else{
    message("Performing median normalization...")
    norm.table <- do.call(rbind, lapply(unique(dataTable$temperature), function(temp){
      # subset dataTable to one temperature
      sub.table <- dataTable[which(dataTable$temperature==temp),]
      norm.list <- sapply(colnames(sub.table),function(coln){
        if (grepl(fcStr, coln)){
          col.median <- median(as.numeric(as.character(sub.table[[coln]])), 
                               na.rm=TRUE)
          norm.col <- as.numeric(as.character(sub.table[[coln]])) / col.median
          return(norm.col)
        }
      })
      norm.df <- data.frame(norm.list[!sapply(norm.list, is.null)])
      
      return(norm.df)
    }))
    colnames(norm.table) <- paste("norm", colnames(norm.table), sep="_")
    big.table <- cbind(dataTable, norm.table)
    message("Done.")
    return(big.table)
  }
}
