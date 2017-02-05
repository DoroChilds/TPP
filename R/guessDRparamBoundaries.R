guessDRparamBoundaries <- function(concentrations){
  ## Compute concentration bounds using largest dilution step among the 
  ## concentration gradient.
  
  # Sort drug concentrations: 
  concentrations <- sort(concentrations, decreasing = FALSE)
  
  # Compute maximal step size:
  conc_diffs   <- diff(concentrations[-1])
  
  dil_step_size <- max(conc_diffs)
  half_step_size <- 0.5*dil_step_size
  
  # Compute boundaries for optimization:
  conc_bds <- c(concentrations[2] - half_step_size, 
                max(concentrations) - half_step_size)
  
  return(conc_bds)
}
