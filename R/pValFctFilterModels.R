pValFctFilterModels <- function(parDF, expNameV, expNameT, minR2, maxPlateau, 
                                flagCheckPlateau){
  ## Check melting curve quality before starting p-value computation by 
  ## applying predefined filters on the curve parameters and R2.
  
  message("Performing quality check on the melting curves of both experiments.")

  r2V <- parDF[, paste("R_sq", expNameV, sep="_")]
  r2T <- parDF[, paste("R_sq", expNameT, sep="_")]
  passedTest1 <- r2V>minR2 & r2T>minR2
  message("1. R2 > ", minR2," (both Experiment): ", sum(passedTest1, na.rm=TRUE)," out of ", length(passedTest1), " models passed.")
  
  if(flagCheckPlateau) {
    ## Only test plateau if the reference group was annotated as "Vehicle" by 
    ## the user:
    plV <- parDF[, paste("plateau", expNameV, sep="_")]
    passedTest2 <- plV<maxPlateau
    message("2. Pl < ", maxPlateau," (Vehicle group only): ", sum(passedTest2, na.rm=TRUE)," out of ", length(passedTest2), " models passed.")    
  } else {
    passedTest2 <- rep(TRUE, length(passedTest1))
  }
            
  passedTest <- passedTest1 & passedTest2
  passedTest[is.na(passedTest)] <- FALSE
  message("=> ", sum(passedTest)," out of ", length(passedTest), " models passed in total and will be used for p-value computation.\n" ) 
  return(passedTest)
}