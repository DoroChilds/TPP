tidyUpReferenceObject <- function(refDatList, refFcColName, refIdColName){
  # Convert reference data from wide to long table:
  refMeasurements <- refDatList$sumResTable$detail %>% tbl_df
  lblsByTemp      <- refDatList$lblsByTemp
  cfg             <- refDatList$tppCfgTable
  
  refTableLong <- convert_trData_wide_to_long(datWide = refMeasurements, 
                                              idColname = refIdColName, 
                                              fcColname = refFcColName,
                                              lblsByTemp= lblsByTemp,
                                              experiments = cfg$Experiment) %>%
    rename_("uniqueID" = refIdColName)
  return(refTableLong)
}