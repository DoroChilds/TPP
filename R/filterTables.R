filterTables <- function(data, normReqs){
  ## Filter all data sets by specified criteria on fold changes and metadata columns

  message("Creating normalization set:")
  
  ## Only regard proteins that occur in all treatment groups:
  grNames <- names(data)
  
  message("\t1. Filtering by non fold change columns:")
  filtersOther <- normReqs[["otherRequirements"]]
  cols <- as.character(filtersOther$colName) # Column to filter
  lb   <- filtersOther$thresholdLower # Lower bounds
  ub   <- filtersOther$thresholdUpper # Upper bounds
  
  if (length(cols) > 0){
    fcListFiltered <- sapply(grNames, function(n){
      message("Filtering by annotation column(s) '", paste(cols, collapse="', '"), "' in treatment group: ", n)
      filterOther(data=data[[n]], cols=cols, lb=lb, ub=ub)},
      simplify=FALSE)
  } else fcListFiltered <- data
  
  ## Determine proteins that occur in each experiment (jointP):
  message("\t2. Find jointP:")
  jointP <- findJointP(data=fcListFiltered)
  fcListJointP   <- sapply(fcListFiltered, function(x){x[jointP,]})
  
  ## Filter by fold changes:
  message("\t3. Filtering fold changes:")
  filtersFC <- normReqs[["fcRequirements"]]
  pos <- filtersFC$fcColumn       # Column positions
  lb  <- filtersFC$thresholdLower # Lower bounds
  ub  <- filtersFC$thresholdUpper # Upper bounds
  
  fcListFiltered <- sapply(grNames, function(n){message("Filtering fold changes in treatment group: ", n)
    filterFCs(data=fcListJointP[[n]], pos=pos, lb=lb, ub=ub)},
    simplify=FALSE)
  
  ## Find treatment group with largest protein set after filtering (normP):
  protNum    <- sapply(fcListFiltered, nrow)
  grNormP    <- grNames[which.max(protNum)]   # Name of treatment group
  namesNormP <- featureNames(fcListFiltered[[grNormP]]) # Protein IDs in normP
  
  message("Experiment with most remaining proteins after filtering: ", grNormP)
  message("-> NormP contains ", length(namesNormP), " proteins.")
  
  return(list("group"       = grNormP,
              "protein_IDs" = namesNormP))
}