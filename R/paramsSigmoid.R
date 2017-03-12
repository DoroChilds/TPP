paramsSigmoid <- function(model, xRange, y){
  ## Compute melting curve and fit parameters for sigmoidal model.
  
  ## 1) Obtain parameter 'a':
  a <- coef(model)[["a"]]
  ## 2) Obtain parameter 'b':
  b <- coef(model)[["b"]]
  ## 3) Obtain plateau parameter:
  pl <- coefficients(model)[["Pl"]]
  ## 4) Fit melting point:
  mp <- meltingPoint(model=model, xRange=xRange)
  ## 5) Find inflection point:
  ip <- inflectionPoint(model=model, xRange=xRange)
  ## 6) Compute slope:
  sl <- meltingCurveSlope(model=model, xInfl=ip)
  ## 7) Compute R2:
  r2 <- rSquared(model=model, y=y)    
  ## 8) Return results:
  res <- c("a"=a, "b"=b, "meltPoint"=mp, "inflPoint"=ip, "slope"=sl, "plateau"=pl, "R_sq"=r2)
  return(res)
}
