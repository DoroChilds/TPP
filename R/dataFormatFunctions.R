eSetsToLongTable_fc <- function(data){
  expNames <- names(data)
  longTabFC <- c()
  for (en in expNames){
    datTmp <- data[[en]]
    
    # Retrieve protein ids, pData, fData and fold changes:
    ids <- featureNames(datTmp)
    pDatWide <- pData(datTmp)
    fcWide  <- data.frame(id=ids, exprs(datTmp))
    
    # Long table of fold changes:
    labelVals <- pDatWide[,2]
    fcLong <- reshape(data=fcWide, direction="long", timevar="labelValue",
                      varying=list(sampleNames(datTmp)), v.names="foldChange",
                      times=labelVals, idvar="id")
    fcLong$labelName <- mapvalues(fcLong$labelValue, from=labelVals, to=pDatWide$label)
    fcLong$experiment <- en
    fcLong$id <- as.character(fcLong$id)
    longTabFC <- rbind(longTabFC, fcLong)
  }
  longTabFC <- arrange(longTabFC, id)
  return(longTabFC)
}

eSetsToLongTable_fData <- function(data){
  expNames <- names(data)
  longTabAnnot <- c()
  for (en in expNames){
    datTmp <- data[[en]]
    
    # Retrieve protein ids, pData, fData and fold changes:
    ids <- featureNames(datTmp)
    fDatWide <- pData(featureData(datTmp))
    
    # Long table of fold changes:
    colnames(fDatWide) <- gsub("([^[:alnum:]])", "_", colnames(fDatWide))
    fDatLong <- reshape(data=data.frame(id=ids, fDatWide), 
                        direction="long", timevar="variable",
                        varying=list(colnames(fDatWide)), v.names="value",
                        times=colnames(fDatWide), idvar="id")
    fDatLong$experiment <- en
    fDatLong$id <- as.character(fDatLong$id)
    longTabAnnot <- rbind(longTabAnnot, fDatLong)
  }
  longTabAnnot <- arrange(longTabAnnot, id)
  return(longTabAnnot)
}


retrieveDataFromESets_TR <- function(data){
  ## Merge row annotation data and fold changes over all experiments
  
  ## 1. Preparation
  expNames<- names(data)
  list1=list2=list4=list3=list5=list6 = vector(mode="list", length=length(expNames))
  names(list1)=names(list2)=names(list4)=names(list3)=names(list5)=names(list6) = expNames
  
  ## 2. Iterate over all experiments and retrieve data
  for (en in expNames){
    setTmp <- data[[en]]
    
    
    ## Split annotation data (stored as featureData in the expressionSets) into
    ## a data frame of melting curve parameters, model information (boolean
    ## variables indicating whether sufficient non-missing values were available
    ## for model fit and whether the model converged successfully).
    fDat        <- pData(featureData(setTmp))
    
    ## Specify column names:
    cols1 <- meltCurveParamNames(returnParNames=TRUE, returnPerformanceInfo=FALSE)
    cols2 <- meltCurveParamNames(returnParNames=FALSE,returnPerformanceInfo=TRUE)
    
    ciOpt <- getOption("TPPTR_CI")
    if(!is.null(ciOpt)){
      if (ciOpt){
        cols2 = c(cols2, "CI_meltPointUpper", "CI_meltPointLower", "CI_meltPoint_delta")
      }
    }    
    cols6 <- "plot"
    cols3 <- setdiff(colnames(fDat), c(cols1,cols2, cols6))
    
    
    ## Split featureData into separate data frames:
    df1  <- as.data.frame(fDat[, cols1])
    df2  <- as.data.frame(fDat[, cols2])
    df3  <- as.data.frame(fDat[, cols3])
    df6  <- subset(fDat, select=cols6)
    
    
    ## Retrieve fold change matrix from current expressionSet and convert to 
    ## data frame:
    df4 <- as.data.frame(exprs(setTmp))
    
    ## Data frame with indicators which proteins where identified in which experiment:
    df5 <- data.frame("protein_identified_in" = rep(TRUE, nrow(setTmp)))
    
    ## Append experiment id to all data frame columns to make them unique when 
    ## combined to big experiment-spanning results table:
    colnames(df1) <- paste(colnames(df1), en, sep="_")
    colnames(df2) <- paste(colnames(df2), en, sep="_")
    colnames(df3) <- paste(colnames(df3), en, sep="_")
    colnames(df4) <- paste(colnames(df4), en, sep="_")
    colnames(df5) <- paste(colnames(df5), en, sep="_")
    colnames(df6) <- paste(colnames(df6), en, sep="_")
    
    ## If data was normalized, add suffix 'norm_' to the fold change column 
    ## names. Normalized data is recognized by the values of the normalization 
    ## coefficients in the fold change column annotation.
    flagIsNormalized <- any(!is.na(pData(data[[en]])$normCoeff))
    if (flagIsNormalized) {
      colnames(df4) <- paste("norm", colnames(df4), sep="_")
    }
    
    ## Add protein ID column so that the data frames of multiple experiment 
    ## (with different subsets of proteins detected in each experiment) can 
    ## later be merged together in a robust way:
    idsTmp <- featureNames(setTmp)
    df1 <- data.frame(Protein_ID=idsTmp, df1, stringsAsFactors=FALSE)
    df2 <- data.frame(Protein_ID=idsTmp, df2, stringsAsFactors=FALSE)
    df3 <- data.frame(Protein_ID=idsTmp, df3, stringsAsFactors=FALSE)
    df4 <- data.frame(Protein_ID=idsTmp, df4, stringsAsFactors=FALSE)
    df5 <- data.frame(Protein_ID=idsTmp, df5, stringsAsFactors=FALSE)
    df6 <- data.frame(Protein_ID=idsTmp, df6, stringsAsFactors=FALSE)
    
    ## Store data frames of each experiment in a list. This will enable
    ## easy and robust merging using plyr::join_all.
    list1[[en]] <- df1
    list2[[en]] <- df2
    list3[[en]] <- df3
    list4[[en]] <- df4
    list5[[en]] <- df5
    list6[[en]] <- df6
  }
  merged1 <- arrange(join_all(list1, by="Protein_ID", type="full"), Protein_ID)
  merged2 <- arrange(join_all(list2, by="Protein_ID", type="full"), Protein_ID)
  merged3 <- arrange(join_all(list3, by="Protein_ID", type="full"), Protein_ID)
  merged4 <- arrange(join_all(list4, by="Protein_ID", type="full"), Protein_ID)
  merged5 <- arrange(join_all(list5, by="Protein_ID", type="full"), Protein_ID)
  merged6 <- arrange(join_all(list6, by="Protein_ID", type="full"), Protein_ID)
  
  ## Insert FALSE if a protein was not present in an experiment (instead of the
  ## NAs generated by the join_all function):
  for (en in expNames){
    name <- paste("protein_identified_in", en, sep="_")
    x    <- merged5[, name]
    x[is.na(x)] <- FALSE
    merged5[, name] <- x
  }
  
  ## Merge plot columns (columns of individual experiments can contain missing 
  ## values if experiment did not provide enough data for plotting):
  plotCols <- grep("plot", colnames(merged6), value = TRUE)
  if (length(plotCols)>0){
    allPlots <- data.frame(Protein_ID = merged6$Protein_ID, "plot"=merge_cols(merged6[,plotCols], collapse='|'))
    merged6 <- join(merged6, allPlots, by="Protein_ID")
  }
  merged6 <- subset(merged6, select = !colnames(merged6) %in% plotCols)
  
  ## Return results:
  return(list(curveParDF   = merged1,
              modelInfoDF  = merged2,
              otherAnnotDF = merged3,
              fcDF         = merged4,
              presenceDF   = merged5,
              plotCol      = merged6))
}


retrieveDataFromESets_CCR <- function(data){
  ## Merge row annotation data and fold changes over all experiments
  
  ## 1. Preparation
  expNames<- names(data)
  list1=list2=list4=list3=list5=list6=list7=list8=list9=list10= 
    vector(mode="list", length=length(expNames))
  names(list1)=names(list2)=names(list4)=names(list3)=names(list5)=names(list6)=
    names(list7)=names(list8)=names(list9)=names(list10) = expNames
  
  ## 2. Iterate over all experiments and retrieve data
  for (en in expNames){
    setTmp <- data[[en]]
    
    ## Split annotation data (stored as featureData in the expressionSets) into
    ## a data frame of curve parameters, model information, fold changes, and
    ## annotation columns.
    fDat <- pData(featureData(setTmp))
    fNames <- colnames(fDat)
    fcNames <- sampleNames(setTmp)
    idsTmp <- featureNames(setTmp)
    
    ## Specify column names:
    cols1 <- drCurveParamNames(names = TRUE, info = FALSE)
    cols2 <- drCurveParamNames(names = FALSE, info = TRUE)
    cols3 <- paste(fcNames, "unmodified", sep="_")
    cols4 <- paste(fcNames, "median_normalized", sep="_")
    cols5 <- paste(fcNames, "transformed", sep="_")
    cols10 <- paste(fcNames, "normalized_to_lowest_conc", sep="_")
    cols6 <- "plot"
    cols9 <- c("compound_effect", "meets_FC_requirement")
    cols7 <- setdiff(fNames, c(cols1,cols2, cols3, cols4, cols5, cols6, cols9, cols10))
    
    ## Split featureData into separate data frames:
    df1  <- subset(fDat, select=cols1)
    df2  <- subset(fDat, select=cols2)
    df3  <- subset(fDat, select=cols3)
    df4  <- subset(fDat, select=cols4)
    df5  <- subset(fDat, select=cols5)
    df6  <- subset(fDat, select=cols6)
    df7  <- subset(fDat, select=cols7)
    df9  <- subset(fDat, select=cols9)
    df10  <- subset(fDat, select=cols10)
    
    ## Data frame with indicators which proteins were identified per experiment:
    df8 <- data.frame("protein_identified_in" = rep(TRUE, nrow(setTmp)))
    
    ## Append experiment id to all data frame columns to make them unique when 
    ## combined to big experiment-spanning results table:
    colnames(df1) <- paste(colnames(df1), en, sep="_")
    colnames(df2) <- paste(colnames(df2), en, sep="_")
    colnames(df3) <- paste(colnames(df3), en, sep="_")
    colnames(df4) <- paste(colnames(df4), en, sep="_")
    colnames(df5) <- paste(colnames(df5), en, sep="_")
    colnames(df6) <- paste(colnames(df6), en, sep="_")
    colnames(df7) <- paste(colnames(df7), en, sep="_")
    colnames(df8) <- paste(colnames(df8), en, sep="_")
    colnames(df9) <- paste(colnames(df9), en, sep="_")
    colnames(df10) <- paste(colnames(df10), en, sep="_")
    
    ## Add protein ID column so that the data frames of multiple experiment 
    ## (with different subsets of proteins detected in each experiment) can 
    ## later be merged together in a robust way:
    df1 <- data.frame(Protein_ID=idsTmp, df1, stringsAsFactors=FALSE)
    df2 <- data.frame(Protein_ID=idsTmp, df2, stringsAsFactors=FALSE)
    df3 <- data.frame(Protein_ID=idsTmp, df3, stringsAsFactors=FALSE)
    df4 <- data.frame(Protein_ID=idsTmp, df4, stringsAsFactors=FALSE)
    df5 <- data.frame(Protein_ID=idsTmp, df5, stringsAsFactors=FALSE)
    df6 <- data.frame(Protein_ID=idsTmp, df6, stringsAsFactors=FALSE)
    df7 <- data.frame(Protein_ID=idsTmp, df7, stringsAsFactors=FALSE)
    df8 <- data.frame(Protein_ID=idsTmp, df8, stringsAsFactors=FALSE)
    df9 <- data.frame(Protein_ID=idsTmp, df9, stringsAsFactors=FALSE)
    df10 <- data.frame(Protein_ID=idsTmp, df10, stringsAsFactors=FALSE)
    
    ## Store data frames of each experiment in a list. This will enable
    ## easy and robust merging using plyr::join_all.
    list1[[en]] <- df1
    list2[[en]] <- df2
    list3[[en]] <- df3
    list4[[en]] <- df4
    list5[[en]] <- df5
    list6[[en]] <- df6
    list7[[en]] <- df7
    list8[[en]] <- df8
    list9[[en]] <- df9
    list10[[en]] <- df10
  }
  merged1 <- arrange(join_all(list1, by="Protein_ID", type="full"), Protein_ID)
  merged2 <- arrange(join_all(list2, by="Protein_ID", type="full"), Protein_ID)
  merged3 <- arrange(join_all(list3, by="Protein_ID", type="full"), Protein_ID)
  merged4 <- arrange(join_all(list4, by="Protein_ID", type="full"), Protein_ID)
  merged5 <- arrange(join_all(list5, by="Protein_ID", type="full"), Protein_ID)
  merged6 <- arrange(join_all(list6, by="Protein_ID", type="full"), Protein_ID)
  merged7 <- arrange(join_all(list7, by="Protein_ID", type="full"), Protein_ID)
  merged8 <- arrange(join_all(list8, by="Protein_ID", type="full"), Protein_ID)
  merged9 <- arrange(join_all(list9, by="Protein_ID", type="full"), Protein_ID)
  merged10 <- arrange(join_all(list10, by="Protein_ID", type="full"), Protein_ID)
  
  ## Insert FALSE if a protein was not present in an experiment (instead of the
  ## NAs generated by the join_all function):
  for (en in expNames){
    name <- paste("protein_identified_in", en, sep="_")
    x    <- merged8[, name]
    x[is.na(x)] <- FALSE
    merged8[, name] <- x
  }
  
  ## Merge plot columns (columns of individual experiments can contain missing 
  ## values if experiment did not provide enough data for plotting):
  plotCols <- grep("plot", colnames(merged6), value = TRUE)
  if (length(plotCols)>0){
    allPlots <- data.frame(Protein_ID = merged6$Protein_ID, "plot"=merge_cols(merged6[,plotCols], collapse='|'))
    merged6 <- join(merged6, allPlots, by="Protein_ID")
  }
  merged6 <- subset(merged6, select = !colnames(merged6) %in% plotCols)
  
  ## Return results:
  return(list(modelPars    = merged1,
              modelInfo    = merged2,
              fcOrig       = merged3,
              fcRefNorm    = merged10,
              fcNorm       = merged4,
              fcTransf     = merged5,
              plotCol      = merged6,
              otherAnnotDF = merged7,
              presenceDF   = merged8,
              transfDF     = merged9))
}

merge_cols <- function(data, fun = paste, ...) {
  # Helper function merge columns of matrix into one vector
  #
  # @param data data matrix
  # @param fun function to deal with different values per row
  # @param ... addtional arguments to fun
  #
  # @return vector with combined columns
     
     data <- as.matrix(data)
     
     res <- apply(data, 1, function(xx) fun(unique(na.omit(xx)), ...))
     
     res[res == ""] <- NA

     res     
}
