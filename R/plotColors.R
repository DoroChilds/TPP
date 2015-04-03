plotColors <- function(expConditions, expReplicates){
  ## Create matching color pairs for vehicle and treatment groups
  groupColors <- FALSE
  if (all(!is.na(expConditions)) && all(!is.na(expReplicates))){
    expReplicates <- as.numeric(expReplicates)
    condLevels <- sort(unique(expConditions))
    repLevels  <- sort(unique(expReplicates))
    condNum <- length(condLevels)
    repNum  <- length(repLevels)
    
    ## Check if several experiments within a replicate got the same condition assigned
    ## (can happen, for example, when all repliacate entries got default value "1")
    if (max(table(paste(expConditions, expReplicates))) == 1){
      ## Check if numbers of conditions and replicates do not exceed maximum
      if (condNum==2 && repNum<=8){ # brewer pal can only produce up to 6 color pairs with 2 intensities each
        if(all.equal(condLevels, c("Treatment", "Vehicle")) & all.equal(repLevels, 1:repNum)){
          groupColors <- TRUE
        }
      }
    }
  }
  
  lightCols <- c("orangered2", "royalblue2", "gray40", "yellow3", "orchid3", "darksalmon") #, "darkseagreen", "yellow3")
  darkCols  <- c("orangered4", "royalblue4", "gray0", "yellow4", "orchid4", "darkorange") #, "lemonchiffon4", "goldenrod4")#, , "blueviolet", )
  if (groupColors==TRUE){
    colorVec <- rep(NA, length(expConditions))
    for (r in repLevels){
      iT <- which(expConditions=="Treatment" & expReplicates==r)
      iV <- which(expConditions=="Vehicle" & expReplicates==r)
      colorVec[iT] <- darkCols[r]
      colorVec[iV] <- lightCols[r]
    }
  } else if (groupColors==FALSE){
    colorVec <- darkCols[1:length(expConditions)]
  }
  return(colorVec)
}