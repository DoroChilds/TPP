#' @title Replace column names for 2D-TPP data
#' 
#' @description Replaces column names for 2D-TPP data so that the TPP-CCR main function can 
#'   deal with the data
#'                            
#' @return A list of dataframe with colnames which match concentrations instead of isobaric 
#'  labels                            
#'                                        
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param data.list list of data frames that contain the data for the 2D-TPP experiment
#' @param intensityStr character string indicating which columns contain the actual 
#'   sumionarea values. Those column names containing the suffix \code{intensityStr} 
#'   will be regarded as containing sumionarea values.
#' @param fcStr character string indicating which columns contain the fold changes
#'   
#' @export
tpp2dReplaceColNames <- function(configTable, data.list, intensityStr, fcStr){
  message("Reformating data...")
  new.list <- lapply(names(data.list), function(l.name){
    # get colnames matching those of the configTable
    coln <- sub(intensityStr, "", colnames(data.list[[l.name]])[grep(intensityStr, 
                colnames(data.list[[l.name]]))])
    # extract Temperature from list name
    temp <- as.numeric(sub(".*_", "", l.name))
    # get matching index from configTable
    ind <- which(configTable$Temperature == temp)
    # get matching concentrations from coln
    col.ids <- which(colnames(configTable) %in% coln)
    concs <- as.character(unlist(configTable[ind,col.ids]))
    intensity.colnames <- paste(intensityStr, concs, sep="")
    colnames(data.list[[l.name]])[grep(intensityStr, colnames(data.list[[l.name]]))] <- intensity.colnames
    # merge fold change columns if they were specified for import
    if (!is.null(fcStr)){
      coln.fc <- sub(fcStr, "", colnames(data.list[[l.name]])[grep(fcStr, colnames(data.list[[l.name]]))]) 
      col.ids.fc <- which(colnames(configTable) %in% coln.fc)
      concs.fc <- as.character(unlist(configTable[ind,col.ids.fc]))
      fc.colnames <- paste(fcStr, concs.fc, sep="")
      colnames(data.list[[l.name]])[grep(fcStr, colnames(data.list[[l.name]]))] <- fc.colnames
    }
    return(data.list[[l.name]])
  })
  #names(new.list) <- names(data.list)
  return(do.call(rbind, new.list))
}