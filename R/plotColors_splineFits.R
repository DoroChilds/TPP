plotColors_splineFits <- function(expConditions){
  # to do: needs unit test
  
  ## Sort condition levels so that a potential 'null model' entry appears
  ## last and gets assigned the black color. This entry was automatically
  ## created by the function 'predict_and_plot_spline_model'
  
  # Check for missing function arguments (manual checks with 'hasArg' would also work)
  expArgs <-c("expConditions")
  fctCall <- match.call()
  myArgs <- names(fctCall)
  sapply(expArgs, function(arg) {
    if (!arg %in% myArgs){
      stop("Error in ", paste(fctCall)[1], 
           ": argument '", arg, "' is missing, with no default", call. = FALSE)
    }
  })
  
  if (any(expConditions == "null model")){
    newLevels <- setdiff(expConditions, "null model")
    addNullColor <- TRUE
  } else {
    newLevels <- expConditions
    addNullColor <- FALSE
  }
  
  if (length(newLevels) == 2){
    newColors <- c("#da7f2d","#808080")
  } else {
    newColors <- plotColors(expConditions = newLevels, 
                            comparisonNums = NA)
  }
  newColors[newLevels == "Vehicle"] = "#808080"
  newColors[newLevels == "Treatment"] = "#da7f2d"
  
  if (addNullColor){
    newLevels <- c(newLevels, "null model")
    newColors <- c(newColors, "black") 
  }
  names(newColors) <- newLevels
  
  return(newColors)
}