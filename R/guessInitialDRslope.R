guessInitialDRslope <- function(dose, response, hill_bds, cpd_effect) {
  ## Guess initial Hill slope for logistic curve fit of TPP-CCR experiments.

  max_slope <- NA
  slopes <- c()
  perc <- 1
  
  # move with sliding window of size 3 over concentrations:
  # * fit line to corresponding responses
  # * calculate slope of fitted line
  for(i in 1:(length(response)-2)) {
    x <- dose[i:(i+2)]
    y <- response[i:(i+2)]
    m = try(fit <- lm(y~x))
    if(class(m) != "try-error"){
      slope <- coef(fit)[2]
      slopes <- c(slopes, slope)
    }
  }
  
  # extract slope:
  # * destabilized proteins: max. slope (>1)
  # * stabilized proteins: min. slope (<-1)
  if(cpd_effect=="stabilized") {
    slopes <- c(slopes[which(slopes<0)], -1)
    max_slope <- quantile(slopes, 1-perc)
  } else {
    slopes <- c(slopes[which(slopes>0)], 1)
    max_slope <- quantile(slopes, perc)
  }
  
  # ensure that guess is within boundaries
  if(max_slope<hill_bds[1]) {
    max_slope <- hill_bds[1]
  } else if(max_slope>hill_bds[2]) {
    max_slope <- hill_bds[2]
  }
  
  
  signif(max_slope,3)
}
