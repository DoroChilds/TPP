fit_spline_model <- function(dat, formula, algorithm){
  if (algorithm == "lm"){
    fitResult <- try(lm(formula, data = dat), silent = TRUE)
  } else if (algorithm == "rlm"){ # Current algorithm for TPP-2D data fitting -> to do: re-use for TPP-2D fits.
    fitResult <- try(rlm(formula, data = dat, maxit = 50), silent = TRUE)
  }
  return(fitResult)
}

