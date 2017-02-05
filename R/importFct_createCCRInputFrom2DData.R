importFct_createCCRInputFrom2DData <- function(configTable, data.list, intensityStr, fcStr){
  ## Replaces column names for 2D-TPP data so that it can be processed by 
  ## analyzeTPPCCR
  ## Returns a list of dataframe with colnames which match concentrations 
  ## instead of isobaric labels
  ## 
  ## @param configTable data frame that specifies important details of the 
  ## 2D-TPP experiment
  ## @param data.list list of data frames that contain the data for the 2D-TPP 
  ## experiment
  ## @param intensityStr character string indicating which columns contain the 
  ##  actual sumionarea values. Those column names containing the suffix 
  ## \code{intensityStr} will be regarded as containing sumionarea values.
  ## @param fcStr character string indicating which columns contain the fold 
  ## changes
  
  message("Reformating data for input into function 'analyzeTPPCCR' ...")
  experiments <- names(data.list)
  
  new.list <- lapply(experiments, function(l.name){
    
    # get colnames matching those of the configTable
    dataset <- data.list[[l.name]]
    allColumns <- colnames(dataset)
    intensityColumns <- grep(intensityStr, allColumns, value = TRUE)
    
    labels <- sub(intensityStr, "", intensityColumns)
    
    # extract Temperature from list name
    temp <- as.numeric(sub(".*_", "", l.name))
    
    # get matching index from configTable
    ind <- which(configTable$Temperature == temp)
    
    # get matching concentrations from labels
    col.ids <- which(colnames(configTable) %in% labels)
    
    concs <- configTable[ind,col.ids] %>% as.numeric %>% as.character
    
    newIntensityColumns <- paste(intensityStr, concs, sep="")
    colnames(dataset)[allColumns %in% intensityColumns] <- newIntensityColumns
    
    # merge fold change columns if they were specified for import
    if (!is.null(fcStr)){
      newColumns <- colnames(dataset)
      coln.fc <- sub(fcStr, "", colnames(dataset)[grep(fcStr, colnames(dataset))]) 
      col.ids.fc <- which(colnames(configTable) %in% coln.fc)
      concs.fc <- configTable[ind,col.ids] %>% as.numeric %>% as.character
      fc.colnames <- paste(fcStr, concs.fc, sep="")
      colnames(dataset)[grep(fcStr, colnames(dataset))] <- fc.colnames
    }
    return(dataset)
  })
  #names(new.list) <- experiments
  return(do.call(rbind, new.list))
}
