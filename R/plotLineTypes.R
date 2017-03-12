plotLineTypes <- function(expConditions){
  ## Create matching linetype pairs for vehicle and treatment groups.
  
  typeVec <- rep(1, length(expConditions))
  ## Assign dashed lines to vehicle curves
  ## From ?aes_linetype_size_shape:
  ## # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 
  ## 6 = twodash
  
  typeVec[expConditions=="Vehicle"] <- 2
  
  return(typeVec)
}