#' @title Calculate fractional abundance and DMSO ratio of successive sumionareas  
#'   (usage of function is only reasonable when at least two temperatures are multiplexed!)
#'   
#' @description Calculates fractional abundance and DMSO ratio of successive sumionareas and 
#'   creates respective columns which are added two the data frame which is handed over
#'  
#' @examples 
#'   data("panobinostat_2DTPP_smallExample")
#'   load(system.file("example_data/2D_example_data/fractAbundInData.RData", package="TPP"))
#'   fracAbund <- tpp2dCalcFractAbundance(configTable = panobinostat_2DTPP_config, 
#'                 dataTable = fractAbundInData,
#'                 intensityStr = "sumionarea_protein_", 
#'                 idVar = "representative")
#'  
#' @param configTable data frame that specifies important details of the TPP-CCR experiment.
#' @param dataTable data frame of TPP-CCR results (e.g. obtained by \code{run2DTPPCCR}).
#' @param intensityStr character string indicating which columns contain the sumionarea values.
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#'  
#' @return Data frame that was handed over with additional columns of fractional abundance
#'   and DMSO1 vs DMSO2 ratio
#'   
#' @export
tpp2dCalcFractAbundance <- function(configTable=NULL, dataTable=NULL, intensityStr=NULL, idVar=NULL){
  message("Calculating fractional abundance and DMSO1 vs. DMSO2 ratio...")
  # determine adjacent temperature pairs
  compound <- as.character(configTable$Compound[1])
  all.temps <- configTable$Temperature
  temps.num <- length(all.temps)
  temp.pairs <- split(all.temps, as.factor(sort(rank(all.temps)%%(temps.num/2))))
  intensity.cols <- grep(intensityStr, colnames(dataTable))
  # loop over all proteins and calculate fractional abundance for each adjacent temperature pair
  subset.list <- lapply(unique(dataTable[[idVar]]), function(prot){
    # subset data according to only one protein
    CCR.subset <- dataTable[which(dataTable[[idVar]] == prot),]
    if (nrow(CCR.subset) == temps.num){
      # loop over all temp pairs and calc DMSO1/DMSO2
      temp.ratios <- sapply(temp.pairs, function(t.pair){
        row1 <- grep(as.character(t.pair[1]), CCR.subset[["temperature"]])
        row2 <- grep(as.character(t.pair[2]), CCR.subset[["temperature"]])
        # calc dmso1 vs dmso2 
        dmso.ratio <- as.numeric(as.character(CCR.subset[row1, intensity.cols[length(intensity.cols)]]))/
          as.numeric(as.character(CCR.subset[row2, intensity.cols[length(intensity.cols)]]))
        # calc sumionarea1 vs sum(sumionareas1, sumionareas2)
        f.abund1 <- 100*sum(as.numeric(sapply(CCR.subset[row1, intensity.cols], as.character)))/
          sum(as.numeric(sapply(CCR.subset[row1, intensity.cols], as.character)), 
              as.numeric(sapply(CCR.subset[row2, intensity.cols], as.character)))
        # calc sumionarea2 vs sum(sumionareas1, sumionareas2)
        f.abund2 <- 100*sum(as.numeric(sapply(CCR.subset[row2, intensity.cols], as.character)))/
          sum(as.numeric(sapply(CCR.subset[row1, intensity.cols], as.character)), 
              as.numeric(sapply(CCR.subset[row2, intensity.cols], as.character)))
        return(list(dmso.ratio, f.abund1, f.abund2))
      })
      CCR.subset$dmso1_vs_dmso2 <- unlist(rep(temp.ratios[1,], each=2))
      CCR.subset$total_sumionarea_fraction <- as.vector(rbind(unlist(temp.ratios[2,]), 
                                                              unlist(temp.ratios[3,])))
    }else{
      CCR.subset$dmso1_vs_dmso2 <- rep(NA, nrow(CCR.subset))
      CCR.subset$total_sumionarea_fraction <- rep(NA, nrow(CCR.subset))
    }
    return(CCR.subset)
  })
  return(do.call(rbind, subset.list))
}