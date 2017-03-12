meltingCurveSlope <- function(model, xInfl){
  ## Compute melting point slope of sigmoidal model.
  if (class(model)!="try-error" & !is.na(xInfl)){
    ## Define expressions for first derivative of sigmoid function:
    strDeriv1  <- fctSigmoidTR(deriv=1)
    exprDeriv1 <- parse(text=strDeriv1)
    
    ## Extract model parameters:
    coeffs <- coefficients(model)
    
    ## Determine slope at inflection point:
    slope <- eval(exprDeriv1, envir=list(a=coeffs[["a"]], b=coeffs[["b"]], 
                                         Pl=coeffs[["Pl"]], x=xInfl))
  } else{slope <- NA}
  return(slope)
}