plot_pVal_distribution <- function(dataWide){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = p_NPARC = p_adj_NPARC = df1 = df2 = df2_moderated = pValType = 
    pValue = label <- NULL
  
  plotDat <- dataWide %>% 
    select(uniqueID, p_NPARC, p_adj_NPARC, df1, df2, df2_moderated) %>%
    gather(pValType, pValue, c(p_NPARC, p_adj_NPARC)) %>%
    mutate(pValType = factor(pValType), 
           df2_moderated = round(df2_moderated, 2)) %>%
    na.omit()
  
  mainTitle = "P-values before and after Benjamini-Hochberg adjustment"
  subTitle = paste("df1 = ", unique(plotDat$df1),
                   ", df2 = ", unique(plotDat$df2),
                   ", df2_moderated = ", unique(plotDat$df2_moderated), sep = "")
  numProt = plotDat %>% group_by(pValType)%>% summarize(n = n()) %>%
    mutate(label = paste("n =", n))
  
  p <- ggplot(data = plotDat, aes(x = pValue)) +
    geom_histogram(aes(y=..density../max(..density..)), # Histogram with density instead of count on y-axis
                   alpha = 1, binwidth = 0.05, na.rm = TRUE) +
    geom_text(data = numProt, aes(x = 0.1, y = 1, label = label), 
              color = "red", inherit.aes = FALSE) +
    ylab("Relative frequency") +
    ylim(c(-0.1,1.1)) + xlim(c(-0.1,1.1)) +
    # geom_density(alpha = .2, fill = "#FF6666") + # Overlay with transparent density plot
    facet_grid(pValType ~ .) +
    ggtitle(paste(mainTitle, subTitle, sep = "\n"))
  return(p)
}

