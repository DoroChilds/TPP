importFct_df_to_eSet <- function(dataframe=NULL, labels, 
                                 labelValues, name, condition, idVar, fcStr, 
                                 qualColName, naStrs, type){
  ## Convert a single TPP-TR or TPP-CCR data frame into an ExpressionSet object.
  message("Importing ",type, " dataset: ", name)
  
  ## -------------------------------------------------------------------------
  ## Make sure that a column with the name stored in idVar exists:
  if (!any(colnames(dataframe) == idVar)){
    stop("idVar Column '", idVar, "' not found in the current dataset.")
  }
  
  ## Replace NAs in ID variable:
  dataframe <- importFct_makeUniqueIDs(inDF=dataframe, idColName=idVar, 
                                       expName=name)
  
  ## Detect fold change column names:
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
  rownames(dataframe) <- dataframe[,idVar]
  
  ## Retrieve fold changes (FC) and sort them by temperature/ concentration:
  fcRaw <- importFct_extractFCs(datDF=dataframe, colsFC=fcCols)
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

importFct_extractFCs <- function(datDF, colsFC){
  datFC           <- subset(datDF, select=colsFC)
  datFCNum        <- colwise(as.numeric)(datFC)
  fcMat           <- as.matrix(datFCNum)
  rownames(fcMat) <- rownames(datFC)
  if (all(is.na(fcMat))){
    stop("Error while loading data: Fold change columns contain only missing values after data import.")
  }
  return(fcMat)
  
}

importFct_fcCols <- function(datDF, fcPrefix, labelSuffix){
  #' Identify fold change columns in the imported data frame.
  #' Fold change columns are composed of the string in 'fcPrefix' (prefix of 
  #' every fold change column) and the strings contained in 'labelSuffix' (the 
  #' TMT labels).
  colsFC    <- paste(fcPrefix, labelSuffix, sep="")
  fcColPos <- match(colsFC, colnames(datDF))
  if (any(is.na(fcColPos))){
    errormsg <- paste("Could not find column(s) '", 
                      paste(colsFC[is.na(fcColPos)], collapse="', '"), 
                      "'. Please speficy the correct fold change column prefix and isotope labels.", 
                      sep="")
    stop(errormsg)
  } else {
    return(colsFC)
  }
}

importFct_removeDuplicates = function(inDF, refColName, nonNAColNames, qualColName=NULL){
  ## Remove duplicate entries in ID column
  message("Removing duplicate identifiers using quality column '", qualColName, "'...")
  nonUniques = unique( inDF[duplicated(inDF[refColName]), refColName] )
  
  retDF = subset(inDF, !(get(refColName) %in% nonUniques))
  
  for(nU in nonUniques){
    tmpDF = subset(inDF, get(refColName) == nU)
    nonNArows = NULL
    for(r in 1:nrow(tmpDF)){
      if(any(!is.na(tmpDF[r, nonNAColNames]))){
        nonNArows = c(nonNArows, r)
      }
    }
    if(length(nonNArows) > 1){
      if(is.null(qualColName)){
        useRow = 1
      } else {
        qualVals = tmpDF[nonNArows, qualColName]
        useRow = match(max(qualVals), qualVals)
      }
    } else {
      useRow = nonNArows[1]
    }
    retDF = rbind(retDF, tmpDF[useRow, ])
  }
  message(nrow(retDF), " out of ", nrow(inDF), " rows kept for further analysis.")
  return(retDF)
}

importFct_makeUniqueIDs = function(inDF, idColName, expName){
  idColumn <- as.character(inDF[,idColName])
  countSuffix <- seq(1:length(idColumn[is.na(idColumn)]))
  idColumn[is.na(idColumn)] <- paste('NA_in', expName, countSuffix, sep='_')
  inDF[,idColName] <- idColumn
  return(inDF)
}

importFct_create_pData <- function(labels, labelValues, fcCols, type){
  if (type == "TR"){
    dataTmp <- data.frame("label"=labels, "temperature"=labelValues, 
                          "normCoeff"=NA, stringsAsFactors=FALSE, 
                          row.names=fcCols)
    metaTmp   <- data.frame(labelDescription = c("Isobaric label",
                                                 "Temperature",
                                                 "Applied normalization coefficient"),
                            row.names = colnames(dataTmp))
  } else if (type == "CCR"){
    dataTmp <- data.frame("label"=labels, "concentration"=labelValues, 
                          "normCoeff"=NA, stringsAsFactors=FALSE, 
                          row.names=fcCols)
    metaTmp   <- data.frame(labelDescription = c("Isobaric label",
                                                 "Concentration",
                                                 "Applied normalization coefficient"),
                            row.names = colnames(dataTmp))
  }
  colInfo <- AnnotatedDataFrame(data=dataTmp, varMetadata=metaTmp)
  return(colInfo)
}

importFct_create_fData <- function(dat, type, fcRaw){
  if (type == "TR"){
    pars <- meltCurveParamNames(TRUE, TRUE)
    
  } else if (type == "CCR"){
    pars <- drCurveParamNames(TRUE, TRUE)
    fcNames <- colnames(fcRaw)
    dat[,paste(fcNames, "unmodified", sep="_")] <- fcRaw
    dat[,paste(fcNames, "normalized_to_lowest_conc", sep="_")] <- NA
    dat[,paste(fcNames, "median_normalized" , sep="_")] <- NA
    dat[,paste(fcNames, "transformed", sep="_")] <- NA
    dat$compound_effect <- NA
    dat$meets_FC_requirement <- NA
  }
  dat[, pars]   <- NA
  dat[, "plot"] <- NA
  
  strMeta <- rep(paste("col", 1:ncol(dat), sep="_"))  
  metaTmp <- data.frame(labelDescription=strMeta, row.names=colnames(dat))
  rowInfo <- AnnotatedDataFrame(data=dat, varMetadata=metaTmp)
  return(rowInfo)
}
