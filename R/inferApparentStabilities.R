inferApparentStabilities <- function(data_2D, 
                                     trRefDataPath, 
                                     idVar = "representative",
                                     fcStr = "norm_rel_fc_protein_",
                                     refFcStr = "norm_rel_fc_protein_"){
  
  message("Normalizing data by TR-reference...")
  
  # Load TPP-TR reference data
  load(trRefDataPath)
  refMeasurements <- tppRefData$sumResTable$detail %>% tbl_df # previous name: detailData
  lblsByTemp <- tppRefData$lblsByTemp
  cfg <- tppRefData$tppCfgTable # tppCfgTable??
  
  # Convert reference data from wide to long table:
  refFcColName <- refFcStr
  refIdColName <- "Protein_ID"
  refTableLong <- convert_trData_wide_to_long(datWide = refMeasurements, 
                                              idColname = refIdColName, 
                                              fcColname = refFcColName,
                                              lblsByTemp= lblsByTemp,
                                              #labels = as.character(lblsByTemp$lbl),
                                              experiments = cfg$Experiment)
  
  # Fit smoothing spline for each protein and evaluate at temperature points
  # of current 2D-TPP data. 
  # This provides normalization coefficients for each protein and temperature.
  newTemperatures <- as.numeric(as.character(unique(data_2D$temperature)))
  model = as.formula("y ~ ns(x, df = 4)")
  normTable <- refTableLong %>% 
    group_by(Protein_ID) %>%
    rename(x = temperature, y = relConc) %>%
    do(fit_and_eval_spline_model(., xNew = newTemperatures, modelFormula = model)) %>%
    rename(temperature = x) %>%
    ungroup
  
  # to do: re-introduce the following checks and error messages:
  # if (length(which(!is.na(protData_reference)))<10){ # check only FC columns instead of all columns?
  #   if (verbose){
  #     message(paste("The TR reference dataset does not supply enough data points for", 
  #                   protID, sep=" ")) 
  #   }
  #   return(NULL)
  # }else if (!is.null(protData_2D) && nrow(protData_2D)>4){
  
  # Convert 2D-TPP data to be analyzed from wide to long table
  tppData_long <- convert_2dData_wide_to_long(datWide = data_2D, 
                                              idColname = idVar,
                                              fcStr = fcStr)
  
  # Normalization to reference:
  tppData_long_normalized <-  tppData_long %>%
    left_join(normTable %>% mutate(Protein_ID = as.character(Protein_ID)), 
              by = c("Protein_ID", "temperature")) %>%
    mutate(fcNormalized = fc * splinePrediction)
  
  message("Done.")
  return(list(normResult = tppData_long_normalized, referenceDataUsed = refTableLong))
}

fit_and_eval_spline_model <- function(df, xNew, modelFormula){
  splineFit <- try(rlm(as.formula(modelFormula), 
                       data = df, maxit = 150),
                   silent = TRUE)
  if (!inherits(splineFit, "try-error")){
    yNew <- predict(splineFit, list(x = xNew))
  } else {
    yNew <- NA
  }
  normDF <- data.frame(x = xNew, splinePrediction = yNew)
  return(normDF)
}

convert_trData_wide_to_long <- function(datWide, idColname, fcColname,
                                        lblsByTemp, experiments){
  # Convert data frame with TPP-TR measurements from wide to long table:
  # Goal:
  # | uniqueID | label | temperature | concentration | relConc| replicate | condition
  datLong <- datWide %>% 
    select_(.dots=c(idColname, grep(fcColname, colnames(.), value = TRUE))) %>% 
    gather_("key", fcColname, grep(fcColname, colnames(.), value = TRUE)) %>% 
    arrange_(idColname) %>%
    mutate(key = gsub(paste(fcColname, "_", sep=""), "", key)) %>%
    rename_(relConc = fcColname)
  
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

convert_2dData_wide_to_long <- function(datWide, idColname, fcStr){
  # Convert 2D-TPP dataset to long table
  datLong <- datWide %>% as.tbl() %>%
    select(matches(paste(idColname, "|temperature|", "^", fcStr, "[0-9,\\.]+$|", "^", fcStr, "[0-9,\\.]+_unmodified", sep = ""))) %>% 
    gather(columnName, fc, contains(fcStr)) %>%
    rename_(Protein_ID = idColname) %>%
    arrange(Protein_ID)
  
  # Add column with drug concentrations
  ptrn <- paste(fcStr, "([0-9,\\.]+)[^0-9]*", sep="")
  oldValues <- unique(datLong$columnName)
  newValues <- paste(sub(ptrn, "\\1", oldValues), "uM", sep="")
  newCol <- plyr::mapvalues(datLong$columnName, oldValues, newValues)
  datLong2 <- datLong %>% mutate(drugConc = newCol)
}