computeNormFactors <- function(data, startPars, maxAttempts, fixedReference=NULL){
  ## Fit sigmoids to median fold changes of each treatment group and determine 
  ## best fit.
  message("Computing normalization coefficients:")
  grNames <- names(data)
  
  ## Compute median curves for each treatment group:
  message("1. Computing fold change medians for proteins in normP.")
  yMat <- sapply(grNames, function(n) {esApply(data[[n]], 2, median, na.rm=TRUE)}, simplify=TRUE)
  
  ## Fit sigmoids to each median curve:
  message("2. Fitting melting curves to medians.")
  xMat <- sapply(grNames, function(n) {data[[n]]$temperature}, simplify=TRUE)
  modelList <- sapply(grNames, function(n){fitSigmoidTR(xVec=xMat[,n], yVec=yMat[,n], startPars=startPars, maxAttempts=maxAttempts)}, simplify=FALSE)
  
  ## Melting curve parametes:
  r2 <- sapply(grNames, function(gn) rSquared(model=modelList[[gn]], y=yMat[,gn]))
  
  if(is.null(fixedReference)){
    ## Select best curve fit:
    refidx <- which.max(r2)
    mBest  <- modelList[[refidx]]
    message(paste("-> Experiment with best model fit: ", grNames[refidx], " (R2: ", signif(r2[refidx],4),")", sep=""))
  } else {
    refidx <- which(grNames == fixedReference)
    mBest <- modelList[[refidx]]
    message(paste("-> Fixed reference experiment: ", grNames[refidx], sep=""))
  }
  
  ## Compute correction factors and store in ExpressionSet object:
  message("3. Computing normalization coefficients")
  y2Mat    <- sapply(grNames, function(n) predict(mBest, newdata=list(x=xMat[,n])))
  dfCoeffs <- y2Mat/yMat
  
  
  ## Return best model:
  return(list("models"      = modelList,
              "medians"     = yMat,
              "tempVals"    = xMat,
              "rSquared"    = r2,
              "bestFit"     = grNames[[refidx]],
              "corrFactors" = dfCoeffs))
}
