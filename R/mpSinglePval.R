mpSinglePval = function(x, r_1, r0, r1){
  ## Compute p-value for a given melting point difference value using 
  ## precomputed distribution quantiles.
  if(!is.na(x)){
    if(x > r0){
      z = (x - r0) / (r1-r0)
      p <- 1/2 * VGAM::erfc(z/sqrt(2))
    } else {
      z = (r0 - x) / (r0-r_1)
      p <- 1/2 * VGAM::erfc(z/sqrt(2))
    }
  } else {
    p <- NA_real_
  }
  return(p)
}
