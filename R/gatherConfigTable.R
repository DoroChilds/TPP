gatherConfigTable <- function(config){
  # to do: needs unit test
  
  ## Convert the user-defined wide config table from a TPP-TR, -CCR or -2D 
  ## experiment into a long table which is easier to use internally.
  ##
  ## Desired column names of the long table:
  ## - experiment: ID of a single multiplexed MS run
  ## - compound: administered compound
  ## - concentration: compound concentration in mol/L
  ## - temperature: administered temperature
  ## - isobaricLabel: isobaric label of a single concentration/temperature pair
  ## - referenceLabel: the isobaric label of the data points used for fold 
  ##                   change calculation
  ## - path: optional path to the data file in 'txt' or 'csv' format. If not 
  ##         empty, there must be exactly one path per experiment.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Experiment = Compound = concentration = Temperature = isobaricLabel = RefCol =
    Path = resultColumn = concStr <- NULL
  
  ## Convert wide config table to long config table:
  allCols <- colnames(config)
  labelCols <- detectLabelColumnsInConfigTable(allColumns = allCols)
  
  cfgLong <- config %>% 
    as.tbl %>%
    gather_("isobaricLabel", "concentration", labelCols)
  
  ## rearrange and rename column names:
  cfgLong <- cfgLong %>% 
    select(experiment = Experiment, 
           compound = Compound,
           concentration, 
           temperature = Temperature,
           isobaricLabel,
           referenceLabel = RefCol,
           path = Path)
  
  ## remove concentrations that were encoded by special characters like '-'.
  ## convert the remaining concentrations to numeric values with unit mol/L:
  cfgLong <- cfgLong %>% 
    filter(grepl("[[:alnum:]|\\.]", concentration)) %>%
    mutate(concentration = as.numeric(concentration) * 1e-6)
  
  ## create and store strings that will later be used to create the column names 
  ## in the wide result table (example: "sumionarea_protein_128H_131L_0uM"):
  colSuffix <- cfgLong %>%
    group_by(concentration) %>%
    summarize(resultColumn = paste(unique(isobaricLabel), collapse = "_")) %>%
    mutate(concStr = paste0(concentration * 1e6, "uM") %>% gsub("-", ".", .)) %>%
    mutate(resultColumn = paste0(resultColumn, "_", concStr))
    
  cfgLong <- cfgLong %>%
    left_join(colSuffix, by = "concentration")
  
  return(cfgLong)
}