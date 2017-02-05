mpSinglePval = function(x, r_1, r0, r1, type){
  ## Compute p-value for a given melting point difference value using 
  ## precomputed distribution quantiles.
  if(!is.na(x)){
    if(x > r0){
      z = (x - r0) / (r1-r0)
      p <- VGAM::erfc(z/sqrt(2))
    } else {
      z = (r0 - x) / (r0-r_1)
      p <- VGAM::erfc(z/sqrt(2))
    }
  } else {
    p <- z <- NA_real_
  }
  if (type == "p"){
    return(p)    
  } else if (type == "z"){
    return(z)
  }
}
