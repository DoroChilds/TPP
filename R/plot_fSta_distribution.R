plot_fSta_distribution <- function(dataLong){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = staType = staValue = df1 = df2 = df2_moderated = DF1 = 
    label <- NULL
  
  plotDat <- dataLong %>% 
    select(uniqueID, staType, staValue, df1, df2, df2_moderated) %>%
    mutate(staType = factor(staType), df2_moderated = round(df2_moderated, 2)) %>%
    na.omit()
  
  mainTitle = "F-statistics and moderated F-statistics"
  subTitle = paste("Type = ", unique(plotDat$staType), sep = "")
  numProt = plotDat %>% group_by(staType, df1, df2, df2_moderated)%>% 
    summarize(n = n(), DF1 = unique(df1), DF2 = unique(df2), MOD.DF2 = unique(df2_moderated)) %>%
    mutate(label = paste("n =", n, ", df1 =", DF1, ", df2 =", df2, ", df2_moderated =", df2_moderated, sep = " "))
  
  p <- ggplot(data = plotDat, aes(x = staValue)) +
    geom_histogram(aes(y=..density../max(..density..)), # Histogram with density instead of count on y-axis
                   alpha = 1, bins = min(nrow(plotDat), 100), na.rm = TRUE) +
    geom_text(data = numProt, aes(x = 0.1, y = 1, label = label), 
              color = "red", inherit.aes = FALSE, hjust = "left", vjust = "top") +
    # geom_density(alpha = .2, fill = "#FF6666") + # Overlay with transparent density plot
    facet_grid(df1 + df2 + df2_moderated ~ .) +
    ylim(c(0,1)) +
    ylab("Relative frequency") + xlab("F statistic") +
    ggtitle(paste(mainTitle, subTitle, sep = "\n"))
  # stat_function(fun = df, colour = "red", args = list(df1 = df1, df2 = df2))
  
  return(p)
}
