importFct_df_to_eSet <- function(dataframe, labels, labelValues, name, 
                                 condition, idVar, fcStr, qualColName, naStrs, 
                                 type){
  ## Convert a single TPP-TR or TPP-CCR data frame into an ExpressionSet object.
  message("\nImporting ",type, " dataset: ", name)
  
  ## -------------------------------------------------------------------------
  ## Make sure that a column with the name stored in idVar exists:
  if (!any(colnames(dataframe) == idVar)){
    stop("idVar Column '", idVar, "' not found in the current dataset.")
  }
  
  ## Replace NAs in ID variable:
  dataframe <- importFct_makeUniqueIDs(inDF=dataframe, idColName=idVar, 
                                       expName=name)
  
  ## Detect fold change column names:
  if (!is.numeric(labelValues)){
    stop("Temperatures or concentrations (provided by argument 'labelValues') must be numeric vectors", call. = TRUE)
  }
  
  o           <- order(labelValues, decreasing=FALSE)
  labelValues <- labelValues[o]
  labels      <- labels[o]
  fcCols <- importFct_fcCols(datDF=dataframe, fcPrefix=fcStr, labelSuffix=labels)
  
  ## Make sure that ID variable 'idvar' is unique:
  dataframe <- importFct_removeDuplicates(inDF=dataframe, 
                                          refColName=idVar, 
                                          nonNAColNames=fcCols, 
                                          qualColName=qualColName)
  
  ## Use ID variable as row names:
  rownames(dataframe) <- dataframe[[idVar]]
  
  ## Retrieve fold changes (FC) and sort them by temperature/ concentration:
  fcRaw <- importFct_extractFCs(datDF=dataframe, colsFC=fcCols, type=type)
  #     fcRefNorm <- importFct_normalizeToReference(foldChanges=fcRaw, 
  #                                                 refCol=which.min(labelValues))
  
  ## Retrieve protein annotation columns:  
  colsAnnot <- setdiff(colnames(dataframe), colnames(fcRaw))
  datAnnot <- subset(dataframe, select=colsAnnot)
  datAnnot[,idVar] <- NULL
  
  ## -------------------------------------------------------------------------
  ## Prepare all objects necessary for ExpressionSet construction
  
  ## 1.) Create AnnotatedDataFrame objects with column annotation (temperature 
  ## conditions for TPP-TR experiment, concentrations for TPP-CCR experiment).
  ## Will be passed to the phenoData slot of the ExpressionSets.
  colInfo <- importFct_create_pData(labels=labels, labelValues=labelValues, 
                                    fcCols=colnames(fcRaw), type=type)
  
  ## 2.) Create AnnotadedDataFrame objects with row annotation (protein IDs and 
  ## further measurements). Will be passed to the featureData slot.
  rowInfo <- importFct_create_fData(dat=datAnnot, type=type, fcRaw=fcRaw)
  
  ## 3.) Create character vector of experiment annotation. Will be stored in the 
  ## annotation slot of the ExpressionSets.
  if (type == "TR"){
    annotStr <- c("name"=name, "condition"=condition)
  } else if (type == "CCR"){
    annotStr <- c("name"=name)
  }
  
  ## -------------------------------------------------------------------------
  ## Construct ExpressionSet:
  fcSet <- ExpressionSet(assayData=fcRaw, phenoData=colInfo, 
                         featureData=rowInfo, annotation=annotStr)
  ## Sort rows:
  fcSet <- fcSet[order(featureNames(fcSet)),]
  
  ## -------------------------------------------------------------------------  
  ## Report size of imported dataset and success rate (how many proteins provide 
  ## enough data to be acutally used for curve fitting?):
  message("  -> ", name, " contains ", nrow(fcSet), " proteins.")
  numValid <- sum(rowSums(!apply(fcRaw, 2, is.na))>2)
  message("  -> ", numValid, " out of ", nrow(fcRaw), " proteins (",
          round(numValid/nrow(fcRaw) * 100,2), 
          "%) suitable for curve fit (criterion: > 2 valid fold changes per protein).")
  
  ## Return result
  return(fcSet)
}

