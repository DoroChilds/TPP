exprSubset <- function(exprSet, subset){
  ## Creates a subset of an ExpressionSet that contains only specified featureNames.
  newSet <- exprSet[featureNames(exprSet) %in% subset, ]
  return(newSet)
}