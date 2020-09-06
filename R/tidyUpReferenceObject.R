tidyUpReferenceObject <- function(refDatList, refFcColName, refIdColName){
  # Convert reference data from wide to long table:
  refMeasurements <- refDatList$sumResTable$detail %>% tibble::as_tibble()
  lblsByTemp      <- refDatList$lblsByTemp
  cfg             <- refDatList$tppCfgTable
  
  refTableLong <- convert_trData_wide_to_long(datWide = refMeasurements, 
                                              idColname = refIdColName, 
                                              fcColname = refFcColName,
                                              lblsByTemp= lblsByTemp,
                                              experiments = cfg$Experiment)
  
  if (refIdColName != "uniqueID"){
    refTableLong[["uniqueID"]] <- refTableLong[[refIdColName]]
    refTableLong[[refIdColName]] <- NULL
  }
  
  return(refTableLong)
}