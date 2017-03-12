inflectionPoint <- function(model, xRange){
  ## Compute inflection point of sigmoidal model.
  if (class(model)!="try-error"){
    ## Define expressions for second derivative of sigmoid function:
    strDeriv2  <- fctSigmoidTR(deriv=2)
    exprDeriv2 <- parse(text=strDeriv2)
    
    ## Extract model parameters:
    coeffs <- coefficients(model)
    
    ## Determine inflection point:
    r <- try(uniroot(f=function(fExpr, a,b,Pl,x){eval(fExpr)},
                     fExpr=exprDeriv2, a=coeffs[["a"]], b=coeffs[["b"]], Pl=coeffs[["Pl"]],
                     interval=xRange, tol=0.0001), silent=TRUE)
    if (class(r) != "try-error"){
      xInfl <- r$root
    } else {xInfl <- NA}
  } else{xInfl <- NA}
  return(xInfl)
}
