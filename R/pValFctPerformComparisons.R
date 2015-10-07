pValFctPerformComparisons <- function(curveParsAllExp, method, controlFilter, 
                                      controlpVal, comparisons){
  ## Preparation:
  testGroups <- comparisons$testGroup
  refGroups  <- comparisons$refGroup
  refIsVehicle <- comparisons$refIsVehicle
  allIDs       <- curveParsAllExp$Protein_ID
  
  ## Start p-values computation for each experiment pair and store results in a 
  ## data frame that can be appended to the final output:
  pValDF <- data.frame(Protein_ID=allIDs, stringsAsFactors=FALSE)
  for (i in 1:length(testGroups)){
    ## Preparation for the current comparison:
    currentRefExp  <- refGroups[i]
    currentTestExp <- testGroups[i]    
    flagVehicleRef <- refIsVehicle[i]
    compName <- paste(currentTestExp, "vs", currentRefExp, sep="_")
    message("Computing p-values for comparison ", compName," ...")
    
    ## Compute melting point differences and minimal slopes for each protein
    compCols <- pValFct_addMPdiff_and_MinSlope_Columns(expNameV=currentRefExp, 
                                                       expNameT=currentTestExp, 
                                                       parDF=curveParsAllExp)
    mpDif <- compCols$mpDiffs
    minSl <- compCols$minSl
    pValDF <- cbind(pValDF, mpDif, minSl)
    
    # Quality filter: R2 > minR2 (both Experiment) + Plateau < maxPl (Vehicle)?
    minr2 <- controlFilter$minR2
    maxpl <- controlFilter$maxPlateau
    passedTest <- pValFctFilterModels(expNameV=currentRefExp, 
                                      expNameT=currentTestExp, 
                                      parDF=curveParsAllExp, minR2=minr2, 
                                      maxPlateau=maxpl,
                                      flagCheckPlateau=flagVehicleRef) 
    idsFiltered  <- allIDs[passedTest]
    mpdiffFiltered <- mpDif[passedTest,1]
    minslFiltered  <- minSl[passedTest,1]
    
    ## Compute p-value for each protein fulfilling the quality check:
    pValsCurrent <- pValFctPerformSingleComparison(minsl=minslFiltered, 
                                                   mpdiff=mpdiffFiltered, 
                                                   method=method, 
                                                   control=controlpVal, 
                                                   comparisonName=compName)
    pValDF_filtered <- data.frame("Protein_ID"=idsFiltered, 
                                  stringsAsFactors=FALSE)
    pValDF_filtered[, paste("pVal_adj", compName, sep="_")] <- pValsCurrent
    
    ## Store current p-values in the output table by matching the filtered IDs 
    ## to the complete set of IDs:
    pValDF <- join(pValDF, pValDF_filtered, by="Protein_ID")
    
    ## Store information about which proteins were used for p-value computation:
    pValDF[, paste("passed_filter", compName, sep="_")] <- passedTest
  }
  return(pValDF)
}