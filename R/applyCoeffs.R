applyCoeffs <- function(data, coeffs){
  ## Normalize fold changes stored in an ExpressionSet object using a given 
  ## vector of coefficients.
  fcOld <- Biobase::exprs(data)
  fcNorm <- t(apply(X=fcOld, MARGIN=1, function(x) x*coeffs))

  Biobase::exprs(data)    <- fcNorm
  data$normCoeff <- coeffs

  return(data)
}