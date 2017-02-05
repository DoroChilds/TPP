convert_2D_cfgTable_to_CCR_cfgTable <- function(configTable){
  ## Create TPP-CCR config files for 2D-TPP experiment
  ##  
  ## @return A config file of type data.frame which can be used for the \code{\link{tpp2dCurveFit}} function
  ##  
  ## @param configTable data frame that specifies important details of the 2D-TPP experiment
  
  if (is.null(configTable) | !is.data.frame(configTable)){
    stop("Please specify a valid configTable of type data.frame")
  }else{
    concentrations <- extractConc(configTable)
    exp.name <- configTable$Compound[1]
    concs <- unique(extractConc(configTable))
    names(concs) <- concs
    out <- data.frame(Experiment = exp.name, 
                      t(concs), 
                      check.names = FALSE, 
                      row.names = c(""))
    return(out) 
  }
}