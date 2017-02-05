plot_2D_data_on_temperature_range <- function(tppData_long_normalized, 
                                              refTableLong){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  Protein_ID = experiment = temperature = relConc = condition = fcNormalized = 
    drugConc <- NULL
  
  message("Plotting splines...")
  # Define color platte for plotting
  nv <- length(unique(tppData_long_normalized$drugConc))
  colors <- colorRampPalette(c("orange", "red", "midnightblue"))(nv)
  
  # Create plot per protein
  allIDs <- unique(tppData_long_normalized$Protein_ID)
  uniqueIDs <- unique(allIDs)
  resultList <- lapply(uniqueIDs, function(id_tmp){
    refDF_tmp <- filter(refTableLong, Protein_ID == id_tmp) %>% 
      mutate(condition = experiment)
    tppDF_tmp <- filter(tppData_long_normalized, Protein_ID == id_tmp)
    
    if (nrow(refDF_tmp) >= 10 && nrow(tppDF_tmp) >= 30 && 
        length(which(is.na(tppDF_tmp$fcNormalized))) < 10){
      tppDF_tmp$drugConc <- as.factor(tppDF_tmp$drugConc)
      tppDF_tmp$drugConc <- factor(tppDF_tmp$drugConc, 
                                   levels=paste0(sort(gsub("[^0-9,\\.]", "",
                                                           levels(tppDF_tmp$drugConc))), "uM"))
      
      p <- ggplot(refDF_tmp, aes(x = temperature, y = relConc)) +
        geom_point(aes(shape=condition), color="gray30") +
        stat_smooth(data = refDF_tmp, method = rlm, method.args = 
                      list(maxit=150),
                    formula=(y ~ ns(x, df=4)), se=FALSE, color="gray30", 
                    linetype=2) +
        geom_point(data = tppDF_tmp, 
                   aes(x = temperature, y = fcNormalized, color = drugConc)) +
        geom_smooth(data = tppDF_tmp, method=lm, 
                    aes(x = temperature, y = fcNormalized, color = drugConc),
                    formula=(y ~ ns(x, df=4)), se=FALSE) +
        scale_colour_manual("concentration", values = colors) +
        xlab("Temperature") +
        ylab("Inferred apparent stability") +
        ggtitle(id_tmp) +
        theme_classic() +
        theme(axis.line.x = element_line(colour = "black", size=0.5, 
                                         linetype="solid"),
              axis.line.y = element_line(colour = "black", size=0.5, 
                                         linetype="solid"))
      return(p)
    } else{
      return(NA) 
    }
  })
  
  message("Done.")
  names(resultList) <- uniqueIDs
  return(resultList)
}