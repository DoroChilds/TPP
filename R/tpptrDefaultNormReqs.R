#' Default filter criteria for fold change normalization
#'
#' Filter criteria as described in the publication.
#'
#' @examples
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=tpptrDefaultNormReqs())
#'
#' @return List with two entries: 'fcRequirements' describes filtering 
#' requirements on fold change columns, 'otherRequirements' contains criteria on 
#' additional metadata columns.
#'
#' @export
tpptrDefaultNormReqs <- function(){
  dfFcFilters <- data.frame("fcColumn"       = c(7  , 9  , 10),
                            "thresholdLower" = c(0.4, 0  , 0),
                            "thresholdUpper" = c(0.6, 0.3, 0.2),
                            stringsAsFactors = FALSE)
  dfOtherFilters <- data.frame(colName        = "qssm", 
                               thresholdLower = 4, 
                               thresholdUpper = Inf,
                               stringsAsFactors = FALSE) 
  filterCrit <- list("fcRequirements"=dfFcFilters, 
                     "otherRequirements"=dfOtherFilters)
  return(filterCrit)
}