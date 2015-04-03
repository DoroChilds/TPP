#' @title Plot minSlope vs. melting point differences
#' @description \code{tpptrQCPlotsMinSlopes_vs_MPdiffs} plots the minSlope vs. 
#'   melting point difference for each protein and highlights proteins with 
#'   significant shifts.
#' @export
#' @return No value returned.
#' @examples
#' data(hdacTR_resultsTable_smallExample)
#' tpptrQCPlotsMinSlopes_vs_MPdiffs(resultTable=resultTable, 
#' expNames=c("Vehicle_1", "Vehicle_2", "Panobinostat_1", "Panobinostat_2"),
#' expRepl=c(1,2,1,2), 
#' expCond=c("Vehicle", "Vehicle", "Treatment", "Treatment"))
#' 
#' @param resultTable Table of the TPP-TR results.
#' @param expNames Names of each experiment.
#' @param ggplotTheme ggplot theme to be used in the plots.
#' @param expRepl Replicate each experiment belongs to.
#' @param expCond Condition (Treatment or Vehicle) each experiment belongs to.
#' @param path Location where to store resulting plot.
tpptrQCPlotsMinSlopes_vs_MPdiffs <- function(resultTable, expNames=NULL, 
                                             ggplotTheme=tppDefaultTheme(),
                                             expRepl=NULL, expCond=NULL, 
                                             path=NULL){
  ## Determine whether comparisons should be made between conditions
  ## (currently only possible for the conditions 'Vehicle' and 'Treatment')
  condLevels <- unique(expCond)
  flagCompareMP <- FALSE
  if (length(condLevels)==2){
    if (all.equal(sort(as.character(condLevels)), c("Treatment", "Vehicle"))){
      flagCompareMP <- TRUE
    }
  }
  
  if (flagCompareMP==TRUE){
    for(r in unique(expRepl)){
      expNameV <- expNames[which(expCond=="Vehicle" & expRepl==r)]
      expNameT <- expNames[which(expCond=="Treatment" & expRepl==r)]
      nameMpDiff    <- paste("diff_meltP", expNameT, "vs", expNameV, sep="_")
      nameMinSlope  <- paste("min_slope", expNameT, "vs", expNameV, sep="_")
      namePVal      <- paste("pVal_adj", expNameT, "vs", expNameV, sep="_")
      
      xMpDiff <- resultTable[,nameMpDiff]
      xMinSl  <- resultTable[,nameMinSlope]
      xpVals  <- resultTable[,namePVal]
      
      plotMinSlopes_vs_MPDiffs(mpDiffs=xMpDiff, minSlopes=xMinSl, pValues=xpVals,
                               plotTheme=ggplotTheme, expName1=expNameT, 
                               expName2=expNameV, path=path)
    }    
    
  } else {
    warning("QC plots of melting point differences can only be created when experiments are assigned to the conditions 'Treatment' and 'Vehicle'.")
  }
  return(NULL)
}