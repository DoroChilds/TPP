#' @title Statistical analysis of melting curve parameters.
#' @description Summarizes the output of a TPP-TR experiment and performs
#'   statistical comparisons of conditions, if appropriate.
#' @details If a TPP-TR experiment was performed with the conditions "Vehicle"
#'   and "Treatment", the melting points between these conditions will be
#'   compared statistically, producing p-values for each protein and replicate.
#'   
#' @return A data frame in which the fit results are stored row-wise for each 
#'   protein.
#'   
#' @examples
#' data(hdacTR_fittedData_smallExample)
#' resultTable <- tpptrResultTable(trDataFitted)
#' subset(resultTable, fulfills_all_4_requirements)$Protein_ID
#' 
#' @param data list of ExpressionSets containing fold changes and metadata. It's
#'   featureData contains the fitted melting curve parameters
#' @param binWidth bin width used for p-value computation
#' @export
#' @references Cox, J., & Mann, M. (2008). MaxQuant enables high peptide 
#'   identification rates, individualized ppb-range mass accuracies and 
#'   proteome-wide protein quantification. Nature biotechnology, 26(12),
#'   1367-1372.

tpptrResultTable <- function(data, binWidth=300){
  message("Creating results table.")
  expNames   <- names(data)
  expConds   <- sapply(data, function(d) d@annotation[["condition"]])
  expRepls   <- as.numeric(sapply(data, function(d) d@annotation[["replicate"]]))
  repLevels  <- unique(expRepls)
  condLevels <- unique(expConds)

  ## Merge row annotation data and fold changes over all experiments:
  listMeltcurvePars <- listModelinfo  <- listFoldchanges <- listOtherAnnotation <- vector(mode="list", length=length(expNames))
  names(listMeltcurvePars) <- names(listModelinfo) <- names(listFoldchanges) <- names(listOtherAnnotation) <- expNames
  for (en in expNames){
    setTmp <- data[[en]]

    ## Split annotation data (stored as featureData in the expressionSets) into
    ## a data frame of melting curve parameters, model information (boolean
    ## variables indicating whether sufficient non-missing values were available
    ## for model fit and whether the model converged successfully).
    dfAnnotDataTmp        <- pData(featureData(setTmp))
    ## Specify column names:
    colnamesCurvePars  <- meltCurveParamNames(returnParNames=TRUE, 
                                              returnPerformanceInfo=FALSE)
    colnamesModelInfo  <- meltCurveParamNames(returnParNames=FALSE, 
                                              returnPerformanceInfo=TRUE)
    colnamesOthers     <- setdiff(colnames(dfAnnotDataTmp), 
                                  c(colnamesCurvePars, colnamesModelInfo, "plot"))
    ## Split featureData into separate data frames:
    dfMeltCurveParsTmp    <- as.data.frame(dfAnnotDataTmp[, colnamesCurvePars])
    dfOtherAnnotationTmp  <- as.data.frame(dfAnnotDataTmp[, colnamesOthers])
    dfModelInfoTmp        <- as.data.frame(dfAnnotDataTmp[, colnamesModelInfo])

    ## Retrieve fold change matrix from current expressionSet and convert to data frame:
    dfFoldchangesTmp           <- as.data.frame(exprs(setTmp))

    ## Append experiment id to all data frame columns to make them unique when combined
    ## to big experiment-spanning results table:
    colnames(dfMeltCurveParsTmp)   <- paste(colnames(dfMeltCurveParsTmp), en, sep="_")
    colnames(dfModelInfoTmp)       <- paste(colnames(dfModelInfoTmp), en, sep="_")
    colnames(dfOtherAnnotationTmp) <- paste(colnames(dfOtherAnnotationTmp), en, sep="_")
    colnames(dfFoldchangesTmp)     <- paste(colnames(dfFoldchangesTmp), en, sep="_")

    ## If data was normalized, add suffix 'norm_' to the fold change column names.
    ## Normalized data is recognized by the values of the normalization coefficients
    ## in the fold change column annotation.
    flagIsNormalized <- any(!is.na(pData(data[[en]])$normCoeff))
    if (flagIsNormalized) {
      colnames(dfFoldchangesTmp) <- paste("norm", colnames(dfFoldchangesTmp), sep="_")
    }
    
    ## Add protein ID column so that the data frames of multiple experiment (with
    ## different subsets of proteins detected in each experiment) can later be
    ## merged together in a robust way:
    protIDsTmp              <- featureNames(setTmp)
    resultColsMeltcurvePars <- data.frame(Protein_ID=protIDsTmp, 
                                          dfMeltCurveParsTmp, 
                                          stringsAsFactors=FALSE)
    resultColsModelInfo     <- data.frame(Protein_ID=protIDsTmp, 
                                          dfModelInfoTmp, 
                                          stringsAsFactors=FALSE)
    resultColsAnnot         <- data.frame(Protein_ID=protIDsTmp, 
                                          dfOtherAnnotationTmp, 
                                          stringsAsFactors=FALSE)
    resultColsFC            <- data.frame(Protein_ID=protIDsTmp, 
                                          dfFoldchangesTmp, 
                                          stringsAsFactors=FALSE)

    ## Store data frames of each experiment in a list. This will enable
    ## easy and robust merging using plyr::join_all.
    listFoldchanges[[en]]     <- resultColsFC
    listMeltcurvePars[[en]]   <- resultColsMeltcurvePars
    listModelinfo[[en]]       <- resultColsModelInfo
    listOtherAnnotation[[en]] <- resultColsAnnot
  }
  mergedFoldchanges   <- join_all(listFoldchanges, by="Protein_ID", type="full")
  mergedMeltcurvePars <- join_all(listMeltcurvePars, by="Protein_ID", type="full")
  mergedModelinfo     <- join_all(listModelinfo, by="Protein_ID", type="full")
  mergedOtherAnnot    <- join_all(listOtherAnnotation, by="Protein_ID", type="full")

  ## Concatenate data for output table:
  outTable <- join(mergedFoldchanges, mergedMeltcurvePars, by="Protein_ID")

  ## Determine whether comparisons should be made between conditions
  ## (currently only possible for the conditions 'Vehicle' and 'Treatment')
  flagCompareMP <- FALSE
  if (length(condLevels)==2){
    if (all.equal(sort(as.character(condLevels)), c("Treatment", "Vehicle"))){
      ## Check if several experiments within a replicate got the same condition assigned
      ## (can happen, for example, when all repliacate entries got default value "1")
      if (max(table(paste(expConds, expRepls))) == 1){
        flagCompareMP <- TRUE        
      } 
    }
  }

  ## Determine whether comparisons should be made between replicates
  ## (currently only possible for exactly 2 replicates)
  flagCompareRep <- FALSE
  if (length(repLevels)==2){
    if (all.equal(repLevels, c(1, 2))){
      flagCompareRep <- TRUE
    }
  }

  ## Compare melting points between Vehicle and Treatment per replicate
  if (flagCompareMP){
    for (r in repLevels){
      expNameV <- expNames[which(expConds=="Vehicle" & expRepls==r)]
      expNameT <- expNames[which(expConds=="Treatment" & expRepls==r)]

      ## Compute melting point differences and minimal slopes for each protein
      mpV <- outTable[, paste("meltPoint", expNameV, sep="_")]
      mpT <- outTable[, paste("meltPoint", expNameT, sep="_")]
      slV <- outTable[, paste("slope", expNameV, sep="_")]
      slT <- outTable[, paste("slope", expNameT, sep="_")]
      mpDif <- mpT - mpV
      minSl <- computeMinimalSlopes(slV=slV, slT=slT)

      ## Store p-values in output table:
      nameMpDiff    <- paste("diff_meltP", expNameT, "vs", expNameV, sep="_")
      nameMinSlope  <- paste("min_slope", expNameT, "vs", expNameV, sep="_")
      statsTmp <- data.frame(Protein_ID=outTable$Protein_ID, 
                             stringsAsFactors=FALSE)
      statsTmp[, nameMpDiff]    <- mpDif
      statsTmp[, nameMinSlope]  <- minSl
      outTable <- join(outTable, statsTmp, by="Protein_ID")

      ## Apply quality filter and compute p-values for each protein fulfilling it:
      r2V <- outTable[, paste("R_sq", expNameV, sep="_")]
      r2T <- outTable[, paste("R_sq", expNameT, sep="_")]
      plV <- outTable[, paste("plateau", expNameV, sep="_")]
      # Test: R2 > 0.8 (both columns) + Plateau < 0.3 (Vehicle)
      passedTest <-resultFilterCurvePars(r2V=r2V, r2T=r2T, plV=plV) 
      idsNew    <- outTable$Protein_ID[passedTest]
      minslNew  <- minSl[passedTest]
      mpdiffNew <- mpDif[passedTest]

      nNew <- length(idsNew)
      if (binWidth > nNew){
        warning(paste("P-value computation for replicate ", r, ": Assigned bin width (",binWidth,") is larger than maximum number of proteins that passed quality control (",nNew,"). Setting it to ", nNew, " instead.", sep=""))
        binWidthNew <- nNew
      } else{
        binWidthNew <- binWidth
      }
      pVals <- computePvalues(ids=idsNew, minSlopes=minslNew, mpDiffs=mpdiffNew, 
                              binWidth=binWidthNew)

      ## Store p-values:
      pvalDF <- data.frame(Protein_ID=idsNew, stringsAsFactors=FALSE)
      namePVal <- paste("pVal_adj", expNameT, "vs", expNameV, sep="_")
      pvalDF[, namePVal] <- pVals
      outTable <- join(outTable, pvalDF, by="Protein_ID")

      ## Store information about which proteins were included in p-value computation:
      outTable[, paste("passed_filter", expNameV, "vs", expNameT, sep="_")] <- passedTest
    }

    if(flagCompareRep){
      ## Summarize Treatment vs. Vehicle comparisons for both replicates
      expNameV1 <- expNames[expConds=="Vehicle" & expRepls==1]
      expNameT1 <- expNames[expConds=="Treatment" & expRepls==1]
      expNameV2 <- expNames[expConds=="Vehicle" & expRepls==2]
      expNameT2 <- expNames[expConds=="Treatment" & expRepls==2]

      ## Melting point difference between both controls
      mpV1 <- outTable[, paste("meltPoint", expNameV1, sep="_")]
      mpV2 <- outTable[, paste("meltPoint", expNameV2, sep="_")]
      mpVehicleDiff <- mpV1 - mpV2
      outTable[, paste("diff_meltP", expNameV1, "vs", expNameV2, sep="_")] <- mpVehicleDiff

      ## Is one of the p values for the two replicate experiments < 0.05 and the
      ## other one < 0.1?
      col1 <- outTable[, paste("pVal_adj", expNameT1, "vs", expNameV1, sep="_")]
      col2 <- outTable[, paste("pVal_adj", expNameT2, "vs", expNameV2, sep="_")]
      # Only sucessful, if both p-values are non-missing values
      pValMin <- pmin(col1, col2, na.rm = FALSE) 
      # Only sucessful, if both p-values are non-missing values
      pValMax <- pmax(col1, col2, na.rm = FALSE) 
      checkCol1 <- pValMin<0.05 & pValMax<0.1

      #outTable[, "min(adj_pVals < 0.05) and max(adj_pVals < 0.1)"] <- checkCol1
      outTable[, "min_pVals_less_0.05_and_max_pVals_less_0.1"] <- checkCol1


      ## Do the melting point shifts in the two control vs treatment experiments
      ## have the same direction (i.e. protein was either stabilized or destabilized
      ## in  both cases)?
      colName1 <- paste("diff_meltP", expNameT1, "vs", expNameV1, sep="_")
      colName2 <- paste("diff_meltP", expNameT2, "vs", expNameV2, sep="_")
      col1 <- outTable[, colName1]
      col2 <- outTable[, colName2]
      checkCol2 <- sign(col1) == sign(col2)
      #newCol <- paste(colName1, "and", colName2, "have the same sign")
      newCol <- "meltP_diffs_have_same_sign"
      outTable[, newCol] <- checkCol2

      ## Are both the melting point differences in the control vs treatment
      ## experiments greater than the melting point difference between the two
      ## untreated controls?
      colT1vsV1 <- paste("diff_meltP", expNameT1, "vs", expNameV1, sep="_")
      colT2vsV2 <- paste("diff_meltP", expNameT2, "vs", expNameV2, sep="_")
      colV1vsV2 <- paste("diff_meltP", expNameV1, "vs", expNameV2, sep="_")

      mpVTDiffs    <- abs(outTable[, c(colT1vsV1, colT2vsV2)])
      mpVTDiffsMin <- apply(mpVTDiffs, 1, min)

      checkCol3 <-  mpVTDiffsMin > abs(mpVehicleDiff)
      #newCol <- paste("min(",colT1vsV1,", ",colT2vsV2,") > ",colV1vsV2, sep="")
      newCol <- "meltP_diffs_T_vs_V_greater_V1_vs_V2"
      outTable[, newCol] <- checkCol3

      ## Is the minimum slope in each of the control vs. treatment experiments < -0.06?
      colT1vsV1 <- paste("min_slope", expNameT1, "vs", expNameV1, sep="_")
      colT2vsV2 <- paste("min_slope", expNameT2, "vs", expNameV2, sep="_")

      minSlopes <- outTable[, c(colT1vsV1, colT2vsV2)]
      checkCol4 <- apply(minSlopes, 1, max) < -0.06
      #newCol <- paste("max(",colT1vsV1,", ", colT2vsV2, ") < -0.06", sep="")
      newCol <- "minSlopes_less_than_0.06"
      outTable[, newCol] <- checkCol4

      ## Does protein fulfill all of the four requirements above?
      checkCol5 <- checkCol1 & checkCol2 & checkCol3 & checkCol4
      newCol <- "fulfills_all_4_requirements" # before: Fulfilling all 4 requirements
      outTable[, newCol] <- checkCol5
    }
  }

  # Add plot column (if applicable)
  flagPlotsExist <- any(!is.na(pData(featureData(data[[1]]))$plot))
  if (flagPlotsExist){
    listPlots <- vector(mode="list", length=length(expNames))
    names(listPlots) <- expNames
    for (en in expNames){
      setTmp <- data[[en]]
      plotsTmp <- featureData(setTmp)$plot
      idsTmp   <- featureNames(setTmp)
      listPlots[[en]] <- data.frame(Protein_ID=idsTmp, Plot=plotsTmp)
    }
    allPlots <- join_all(listPlots, by=c("Protein_ID", "Plot"))
    outTable <- join(outTable, allPlots, by="Protein_ID")
  }

  ## Add columns with additional quality information about model fit (boolean
  ## variables indicating whether sufficient non-missing values were
  ## available for model fit and whether the model converged successfully)
  outTable <- join(outTable, mergedModelinfo, by="Protein_ID")

  ## Add futher annotation columns that were directly imported from the input files:
  outTable <- join(outTable, mergedOtherAnnot, by="Protein_ID")

  ## Sort result table alphabetically according to protein ids:
  outTable <- arrange(outTable, outTable$Protein_ID)
  message("Results table created successfully!\n")
  return(outTable)
}