## ------------------------------------------------------------------------- ##
## function 'applyCoeffs':
## ------------------------------------------------------------------------- ##
test_that(desc="testApplyCoeffs", code={
  x <- matrix(rnorm(10), nrow=2, ncol=5)
  eSet <- Biobase::ExpressionSet(assayData=x)

  normCoeffs <- rnorm(5)
  eSetNormed <- TPP:::applyCoeffs(eSet, normCoeffs)
  
  xNew <- unname(Biobase::exprs(eSetNormed))
  xRef <- rbind(x[1,] * normCoeffs, x[2,] * normCoeffs)
})

