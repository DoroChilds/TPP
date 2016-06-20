guessDRparamBoundaries <- function(concentrations){
  ## Compute concentration bounds using largest dilution step among the concentration gradient.

  conc_diffs   <- diff(concentrations[-1])
  dil_step_size <- max(conc_diffs)
  conc_bds <- c(concentrations[2]-0.5*dil_step_size, concentrations[length(concentrations)]-0.5*dil_step_size)
  return(conc_bds)
}