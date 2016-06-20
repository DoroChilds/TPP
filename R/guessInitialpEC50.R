guessInitialpEC50 <- function(dose, response, conc_bds) {
  ## Guess initial pEC50 value for logistic curve fit of TPP-CCR experiments.
  ## weight concentration by proximity of corresponding response to 0.5
  
  w <- 1/abs(response - 0.5)^100
  w[which(w==Inf)] <- 10^300 # set weight to very high value if response = 0.5
  pec50_init <- signif(sum(w * dose, na.rm=TRUE)/sum(w, na.rm=TRUE),3)
  
  # set pEC50 guess to be within concentration range
  if(pec50_init < conc_bds[1]) {
    pec50_init <- conc_bds[1]
  }
  if(pec50_init > conc_bds[2]) {
    pec50_init <- conc_bds[2]
  }
  
  pec50_init
}