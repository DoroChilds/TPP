inferApparentStabilities <- function(data_2D, dataRef, refIDVar, refFcStr){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  tppRefData = uniqueID = temperature = relConc = x = fc = 
    splinePrediction <- NULL
  
  message("Normalizing data by TR-reference...")
  
  # Obtain settings used for data import (stored as attribute of imported data):
  importSettings <- attr(data_2D, "importSettings")
  idVar <- checkAndReturnDataSetting(importSettings, "proteinIdCol", colnames(data_2D))
  
  # Choose correct fold change column prefix (automatically detects whether
  # to use the prefix for normalized columns).
  finalFcPrefix <- obtain_fcStr_from_df_annotation(dat = data_2D)
  
  # Load TPP-TR reference data
  if (is.character(dataRef)){
    if (file.exists(dataRef)){
      load(dataRef)
      dataRef <- tppRefData
    } else {
      stop("Reference data file ", dataRef, " could not be found.")
    }
  }
  ## Check if reference data is a list produced by the function 
  ## tpp2dTRReferenceObject, or already in a tidy format
  if (is.list(dataRef)){
    dataRef <- tidyUpReferenceObject(refDatList = dataRef, 
                                     refFcColName = refFcStr, 
                                     refIdColName = refIDVar)
  }
  

  # Fit smoothing spline for each protein and evaluate at temperature points
  # of current 2D-TPP data. 
  # This provides normalization coefficients for each protein and temperature.
  newTemperatures <- as.numeric(as.character(unique(data_2D$temperature)))
  model = as.formula("y ~ ns(x, df = 4)")
  normTable <- dataRef %>% 
    group_by(uniqueID) %>%
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
                                              fcStr = finalFcPrefix)
  
  # Normalization to reference:
  tppData_long_normalized <-  tppData_long %>%
    left_join(normTable %>% mutate(uniqueID = as.character(uniqueID)), 
              by = c("uniqueID", "temperature")) %>%
    mutate(fcNormalized = fc * splinePrediction)
  
  message("Done.")
  return(list(normResult = tppData_long_normalized, referenceDataUsed = dataRef))
}

