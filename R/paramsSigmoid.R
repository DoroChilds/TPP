paramsSigmoid <- function(model, xRange, y){
  ## Compute melting curve and fit parameters for sigmoidal model
  ## 1) Fit melting point:
  mp <- meltingPoint(model=model, xRange=xRange)
  ## 2) Find inflection point:
  ip <- inflectionPoint(model=model, xRange=xRange)
  ## 3) Compute slope:
  sl <- meltingCurveSlope(model=model, xInfl=ip)
  ## 4) Obtain plateau parameter:
  pl <- coefficients(model)[["Pl"]]
  ## 5) Compute R2:
  r2 <- rSquared(model=model, y=y)    
  ## 6) Return results:
  res <- c("meltPoint"=mp, "inflPoint"=ip, "slope"=sl, "plateau"=pl, "R_sq"=r2)
  return(res)
}