## ------------------------------------------------------------------------- ##
## function 'applyCoeffs':
## ------------------------------------------------------------------------- ##
test_that(desc="testApplyCoeffs", code={
  cat("Checking function 'applyCoeffs': ")
  x <- matrix(rnorm(10), nrow=2, ncol=5)
  eSet <- ExpressionSet(assayData=x)

  normCoeffs <- rnorm(5)
  eSetNormed <- applyCoeffs(eSet, normCoeffs)
  
  xNew <- unname(exprs(eSetNormed))
  xRef <- rbind(x[1,] * normCoeffs, x[2,] * normCoeffs)
  cat(paste("   ", expect_equal(xNew, xRef)))
})