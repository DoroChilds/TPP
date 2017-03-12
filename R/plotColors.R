plotColors <- function(expConditions, comparisonNums){
  ## Create matching color pairs for vehicle and treatment groups.
  groupColors <- FALSE
  
  if (all(!is.na(expConditions))){
    condLevels <- sort(unique(expConditions))
    compLevels  <- sort(unique(comparisonNums))
    condNum <- length(condLevels)
    compNum  <- length(compLevels)
    
    ## Check if several experiments within a comparison got the same condition assigned
    ## (can happen, for example, when comparing two vehicle experiments against each other)
    if (max(table(paste(expConditions, comparisonNums))) == 1){
      ## Check if numbers of conditions and comparisons do not exceed maximum
      if (condNum==2 && compNum<=8  && compNum>0){ # brewer pal can only produce up to 8 color pairs with 2 intensities each
        if(all.equal(condLevels, c("Treatment", "Vehicle")) & all.equal(compLevels, 1:compNum)){
          groupColors <- TRUE
        }
      }
    }
  }
  
  plotClrsLight <- brewer.pal(n=8, name = "Set2")
  plotClrsDark  <- brewer.pal(n=8, name = "Dark2")
  
  NAsInComparisonNums <- is.na(comparisonNums)
  numNAsInComparisonNums <- sum(NAsInComparisonNums)
  
  if (groupColors==TRUE){
    lightCols = colorRampPalette(plotClrsLight)(length(expConditions) + numNAsInComparisonNums)
    darkCols = colorRampPalette(plotClrsDark)(length(expConditions) + numNAsInComparisonNums)
    
    colorVec <- rep(NA, length(expConditions))
    for (r in compLevels){
      iT <- which(expConditions=="Treatment" & comparisonNums==r)
      iV <- which(expConditions=="Vehicle" & comparisonNums==r)
      colorVec[iT] <- darkCols[r]
      colorVec[iV] <- lightCols[r]
    }
    
    if(numNAsInComparisonNums > 0){
      colorVec[is.na(colorVec)] = darkCols[sum(!is.na(colorVec)) : length(darkCols)]
    }
    
  } else if (groupColors==FALSE){
    colorVec <- colorRampPalette(brewer.pal(n=8, name="Dark2"))(length(expConditions))
  }
  return(colorVec)
}