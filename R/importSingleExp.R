importSingleExp <- function(dataframe=NULL, filename=NULL, labels, labelValues,
                            name, condition, replicate, idVar, 
                            fcStr, qualColName, naStrs, type){
  ## Import a single TPP-TR or -CCR experiment from file or data frame and convert 
  #' it into an ExpressionSet.
  ## -------------------------------------------------------------------------
  ## Obtain measurements from file or data frame
  message("Importing ",type, " dataset: ", name)
  
  ## Check that data is specified either as dataframe or filename, but not both
  isDF <- !is.null(dataframe)
  isF  <- !is.null(filename)
  isBoth <- isDF & isF
  isNone <- !(isDF | isF)
  if (isBoth) stop("Data import function received a filename AND a dataframe object. Please specify only one.")
  if (isNone) stop("Data import function requires a filename or a dataframe object. Please specify one.")
  
  ## If check was sucessful, continue with data import
  if (isF){
    ## If data source is specified as filename, import it into R:
    dataframe <- read.delim(filename, as.is=TRUE, na.strings=naStrs)
  }
  
  ## -------------------------------------------------------------------------
  ## Extract all necessary information from the data table
  
  ## 1.) Determine fold change (FC) and non-FC columns:
  colsFC    <- paste(fcStr, labels, sep="")
  fcColPos <- match(colsFC, colnames(dataframe))
  if (any(is.na(fcColPos)) || (length(fcColPos) != length(labels))){
    stop("Looking for fold change colunms that start with prefix '", fcStr,"', followed by the specified isotope labels.\nUnfortunately, such columns were not found in the data table.\nPlease speficy the correct prefix and isotope labels when starting the package functions 'analyzeTPPTR'/'analyzeTPPCCR' (to run the whole analysis) or 'tpptrImport'/'tppccrImport' (to simply import the data and start further analyses later).")
  }
  colsAnnot <- setdiff(colnames(dataframe), colsFC)
  
  ## 2.) Preprocess table
  ## Replace NAs in ID variable:
  dataframe[,idVar] <- as.character(dataframe[,idVar])
  idColumn <- dataframe[,idVar]
  countSuffix <- seq(1:length(idColumn[is.na(idColumn)]))
  idColumn[is.na(idColumn)] <- paste('NA_in', name, countSuffix, sep='_')
  dataframe[,idVar] <- idColumn
  
  ## Make sure that ID variable 'idvar' is unique:
  dataframe <- removeDuplicates(inDF=dataframe, refColName=idVar, nonNAColNames=colsFC, qualColName=qualColName)
  
  ## Use ID variable as row names:
  idsReduced <- dataframe[,idVar]
  rownames(dataframe)  <- idsReduced
  
  ## 3.) Extract annotation columns:
  datAnnot <- subset(dataframe, select=colsAnnot)
  datAnnot[,idVar] <- NULL
  
  ## 4.) Extract fold changes:
  datFC <- subset(dataframe, select=colsFC)
  datFCNum <- colwise(as.numeric)(datFC)
  fcMat <- as.matrix(datFCNum)
  rownames(fcMat) <- rownames(datFC)
  if (all(is.na(fcMat))){
    stop("Error while loading data: Fold change columns contain only missing values after data import.")
  }
  
  ## -------------------------------------------------------------------------
  ## Prepare all objects necessary for ExpressionSet construction
  
  ## 1.) Create AnnotadedDataFrame objects with column annotation(temperature 
  ## conditions for TPP-TR experiment, concentrations for TPP-CCR experiment).
  ## Will be passed to the phenoData slot of the ExpressionSets. 
  labelDFtmp <- data.frame("label"=labels, "temperature"=labelValues, "normCoeff"=NA,
                           stringsAsFactors=FALSE, row.names=colnames(fcMat))
  metaData   <- data.frame(labelDescription = c("Isotope label",
                                                "Temperature (TPP-TR experiment) or concentration (TPP-CCR experiment)",
                                                "Applied normalization coefficient"),
                           row.names = c("label", "temperature", "normCoeff"))
  if (type == "CCR"){
    colnames(labelDFtmp) <- gsub("temperature", "concentration", colnames(labelDFtmp))
    rownames(metaData) <- gsub("temperature", "concentration", rownames(metaData))
  }
  colInfo    <- AnnotatedDataFrame(data=labelDFtmp, varMetadata=metaData)
  
  ## 2.) Create AnnotadedDataFrame objects with row annotation (protein IDs and further 
  ## measurements). Will be passed to the featureData slot of the ExpressionSets. 
  if (type == "TR"){
    meltCurveParNames <- c(meltCurveParamNames(), "plot")
    placeholderMeltCurvePars <- data.frame(matrix(nrow=nrow(dataframe), ncol=length(meltCurveParNames),
                                                  dimnames=list(idsReduced, meltCurveParNames)))
    datAnnotWithMCPars <- cbind(placeholderMeltCurvePars, datAnnot)
    strMeta   <- c(rep("Melting curve parameter", ncol(placeholderMeltCurvePars)),
                   rep("Column from input data", ncol(datAnnot)))
    metaData  <- data.frame(labelDescription=strMeta, row.names=colnames(datAnnotWithMCPars))
    rowInfo   <- AnnotatedDataFrame(data=datAnnotWithMCPars, varMetadata=metaData)      
  } else if (type == "CCR"){
    datAnnotExtended <- datAnnot
    datAnnotExtended$CompoundEffect <- NA
    strMeta   <- c(rep("Column from input data", ncol(datAnnot)),
                   "Classification of compound effect (stabilizing or destabilizing)")
    metaData  <- data.frame(labelDescription=strMeta, row.names=colnames(datAnnotExtended))
    rowInfo   <- AnnotatedDataFrame(data=datAnnotExtended, varMetadata=metaData)
  }
  
  ## 3.) Create character vector of experiment annotation. Will be stored in the 
  ## annotation slot of the ExpressionSets.
  if (type == "TR"){
    ## Assign default values in case the corresponding columns were not 
    ## specified in the config table:
    if (is.null(condition)) condition <- NA_character_ 
    if (is.null(replicate)) replicate <- NA_character_
    ## Create annotation vector:
    annotStr <- c("name"=name, "condition"=condition, "replicate"=replicate)
  } else if (type == "CCR"){
    ## Create annotation vector:
    annotStr <- c("name"=name)
  }
  
  ## -------------------------------------------------------------------------
  
  ## Construct ExpressionSet:
  fcSet <- ExpressionSet(assayData=fcMat, phenoData=colInfo, featureData=rowInfo, annotation=annotStr)
  ## Sort rows according to protein ID:
  fcSet <- fcSet[order(featureNames(fcSet)),]
  if (type == "TR"){
    ## Sort fold change columns according to temperature (TPP-TR):
    fcSet <- fcSet[,order(fcSet$temperature)]     
  } else if (type == "CCR"){
    ## Sort fold change columns according to concentrations (TPP-CCR):
    fcSet <- fcSet[,order(fcSet$concentration)]
  }
  
  ## -------------------------------------------------------------------------
  
  ## Report size of imported dataset and success rate (how many proteins provide 
  ## enough data to be acutally used for curve fitting?):
  message("  -> ", name, " contains ", nrow(fcSet), " proteins.")
  numValid <- sum(rowSums(!apply(fcMat, 2, is.na))>2)
  message("  -> ", numValid, " out of ", nrow(fcMat), " proteins (",round(numValid/nrow(fcMat) * 100,2), "%) suitable for curve fit (criterium: > 2 valid fold changes per protein).")
  
  ## Return result
  return(fcSet)
}
