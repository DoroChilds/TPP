importCheckExperimentNames <- function(expNames, dataframes){
  ## If data is given as a list of dataframes, check whether the names are
  ## consistent with the 'Experiment' column in configTable (argument expNames):
  if (is.list(dataframes) && !is.data.frame(dataframes)){
    if (!identical(expNames, names(dataframes))){
      stop("The names of the data objects in 'data' differ from the names given 
           in the Experiment column of 'configTable'.")
    }
  }
}