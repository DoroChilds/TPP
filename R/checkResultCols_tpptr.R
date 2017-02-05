checkResultCols_tpptr <- function(pValDF, curveParDF, comparisons){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  refIsVehicle = testIsVehicle <- NULL
  
  ## Currently only applicable for exactly 2 treatments and 2 controls
  
  ## Ignore Vehicle vs. Vehicle comparisons:
  comparisons <- subset(comparisons, !(refIsVehicle & testIsVehicle))
  
  ## Test if number of experiments equals two and if the experiments were 
  ## defined as treatment and vehicle groups:
  flagExpNum <- nrow(comparisons)==2
  flagExpConds <- all(c(comparisons$refIsVehicle, comparisons$testIsTreatment))
  
  outTable <- data.frame(Protein_ID=pValDF$Protein_ID, stringsAsFactors=FALSE)
  if(flagExpNum & flagExpConds){
    ## Summarize both Treatment vs. Vehicle comparisons
    expNameV1 <- comparisons$refGroup[1]
    expNameV2 <- comparisons$refGroup[2]
    expNameT1 <- comparisons$testGroup[1]
    expNameT2 <- comparisons$testGroup[2]
    
    ## Melting point difference between both controls
    mpV1 <- curveParDF[, paste("meltPoint", expNameV1, sep="_")]
    mpV2 <- curveParDF[, paste("meltPoint", expNameV2, sep="_")]
    mpVehicleDiff <- mpV1 - mpV2
    outTable[, paste("diff_meltP", expNameV1, "vs", expNameV2, sep="_")] <- mpVehicleDiff
    
    ## ---------------------------------------------------------------------- ##
    ## 1. Check: Is one of the p values < 0.05 and the other one < 0.1?
    col1 <- pValDF[, paste("pVal_adj", expNameT1, "vs", expNameV1, sep="_")]
    col2 <- pValDF[, paste("pVal_adj", expNameT2, "vs", expNameV2, sep="_")]
    # Only sucessful, if both p-values are non-missing values
    pValMin <- pmin(col1, col2, na.rm = FALSE) 
    # Only sucessful, if both p-values are non-missing values
    pValMax <- pmax(col1, col2, na.rm = FALSE) 
    checkCol1 <- pValMin<0.1 & pValMax<0.2
    outTable[, "min_pVals_less_0.1_and_max_pVals_less_0.2"] <- checkCol1
    
    ## ---------------------------------------------------------------------- ##
    ## 2. Check: Do the melting point shifts in the two control vs treatment 
    ## experiments have the same direction (i.e. protein was either stabilized 
    ## or destabilized in  both cases)?
    col1 <- pValDF[, paste("diff_meltP", expNameT1, "vs", expNameV1, sep="_")]
    col2 <- pValDF[, paste("diff_meltP", expNameT2, "vs", expNameV2, sep="_")]
    checkCol2 <- sign(col1) == sign(col2)
    outTable[, "meltP_diffs_have_same_sign"] <- checkCol2
    
    ## ---------------------------------------------------------------------- ##
    ## 3. Check: Are both the melting point differences in the control vs 
    ## treatment experiments greater than the melting point difference between 
    ## the two untreated controls?
    col1 <- pValDF[, paste("diff_meltP", expNameT1, "vs", expNameV1, sep="_")]
    col2 <- pValDF[, paste("diff_meltP", expNameT2, "vs", expNameV2, sep="_")]
    minAbsDiff <- pmin(abs(col1), abs(col2), na.rm=FALSE)
    checkCol3  <- minAbsDiff > abs(mpVehicleDiff)
    outTable[, "meltP_diffs_T_vs_V_greater_V1_vs_V2"] <- checkCol3
    
    ## ---------------------------------------------------------------------- ##
    ## 4. Check: Is the minimum slope in each of the control vs. treatment 
    ## experiments < -0.06?
    col1 <- pValDF[, paste("min_slope", expNameT1, "vs", expNameV1, sep="_")]
    col2 <- pValDF[, paste("min_slope", expNameT2, "vs", expNameV2, sep="_")]
    maxVal    <- pmax(col1, col2, na.rm=FALSE)
    checkCol4 <- maxVal < -0.06
    outTable[, "minSlopes_less_than_0.06"] <- checkCol4
    
    ## ---------------------------------------------------------------------- ##
    ## 5. Check: Does protein fulfill all of the four requirements above?
    checkCol5 <- checkCol1 & checkCol2 & checkCol3 & checkCol4
    outTable[, "fulfills_all_4_requirements"] <- checkCol5
  }
  return(outTable)
}

