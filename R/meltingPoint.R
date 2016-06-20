meltingPoint <- function(model, xRange){
  ## Compute melting point of sigmoidal model.
  
  if (class(model)=="try-error"){
    return(NA)
  } else{
    meltPStr  <- paste(fctSigmoidTR(deriv=0), "-0.5")
    meltPExpr <- parse(text=meltPStr)
    coeffs <- coefficients(model)
    r <- try(uniroot(f=function(fExpr, a,b,Pl,x){eval(fExpr)},
                     fExpr=meltPExpr, a=coeffs[["a"]], b=coeffs[["b"]], Pl=coeffs[["Pl"]],
                     interval=xRange, tol=0.0001), silent=TRUE)
    
    if (class(r) != "try-error"){
      return(r$root)
    } else {
      return(NA)
    }
  }
}
