generateColors4Temps <- function(configTable, paletteName){
  temp.len <- length(unique(configTable$Temperature))
  color.vals <- rep(NA, temp.len)
  if((temp.len %% 2)==0){
    color.vals[1:(temp.len/2)] <- rev(brewer.pal(n=11, name=paletteName))[1:(temp.len/2)]
    color.vals[((temp.len/2)+1):temp.len] <- 
      rev(brewer.pal(n=11, name=paletteName))[(((temp.len/2)+1) -(temp.len-11)):11]
  }else{
    color.vals[1:round((temp.len/2)+1)] <- rev(brewer.pal(n=11, name=paletteName))[1:round((temp.len/2))]
    #color.vals[round((temp.len/2))+1] <- "#B3B3B3"
    color.vals[round((temp.len/2))+2:temp.len] <- 
      rev(brewer.pal(n=11, name=paletteName))[((round(temp.len/2)+1)+2 - (temp.len-11) ):11]
  }
  names(color.vals) <- unique(configTable$Temperature)
  return(color.vals)
}