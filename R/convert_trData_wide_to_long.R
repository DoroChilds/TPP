convert_trData_wide_to_long <- function(datWide, idColname, fcColname,
                                        lblsByTemp, experiments){
  # Convert data frame with TPP-TR measurements from wide to long table.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  key = lbl = temp = temperature = label = experiment <- NULL
  
  # Goal:
  # | uniqueID | label | temperature | concentration | relConc| replicate | condition
  datLong <- datWide %>% 
    select_(.dots=c(idColname, grep(fcColname, colnames(.), value = TRUE))) %>% 
    gather_("key", fcColname, grep(fcColname, colnames(.), value = TRUE)) %>% 
    arrange_(idColname) %>%
    mutate(key = gsub(paste(fcColname, "_", sep=""), "", key))
  
  if (fcColname != "relConc"){
    datLong[["relConc"]] <- datLong[[fcColname]]
    datLong[[fcColname]] <- NULL
  }
  
  
  # Define separate columns for label and experiment using regular expression,
  # and add column with temperatures per label:
  ptrn_labelNames <- paste(as.character(lblsByTemp$lbl), collapse = "|")
  ptrn_expNames <- paste(experiments, collapse = "|")
  ptrn_final <- paste("(", ptrn_labelNames, ")_(", ptrn_expNames,")", sep = "")
  
  labelInfoTable <- lblsByTemp %>% 
    dplyr::rename(label = lbl, temperature = temp) %>%
    mutate(temperature = as.numeric(as.character(temperature)))
  
  datLong2 <- datLong %>% 
    extract(key, c("label", "experiment"), ptrn_final) %>%
    mutate(label = as.factor(label), experiment = as.factor(experiment)) %>%
    left_join(labelInfoTable, by = "label")
  
  return(datLong2)
}