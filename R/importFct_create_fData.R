importFct_create_fData <- function(dat, type, fcRaw){
  ## Create feature data for the expression sets.
  
  if (type == "TR"){
    pars <- meltCurveParamNames(returnParNames = TRUE, 
                                returnPerformanceInfo = TRUE)
    
  } else if (type == "CCR"){
    pars <- drCurveParamNames(names = TRUE, info = TRUE)
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
