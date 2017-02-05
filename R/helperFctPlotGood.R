helperFctPlotGood <- function(configTable, dataTable, idVar, fcStr,  verbose, 
                              paletteName){
  
  # pre-define global variables
  variable <- NULL
  value <- NULL
  temperature <- NULL
  sig.vals <- NULL
  info <- NULL
  
  compound <- configTable$Compound %>% unique
  ## remove the compound name from the result table columns. To avoid confusion, 
  ## throw an error if more than one compound names were in use # new :
  if (length(compound) == 1){
    oldCols <- colnames(dataTable)
    patternToRemove <- paste0("_", gsub("[^[:alnum:]]", "_", compound))
    newCols <- gsub(patternToRemove, "", oldCols)
    colnames(dataTable) <- newCols
  } else {
    stop("Can currently only handle config tables with one compound name.
         Compound names detected in your config table: '", 
         paste(compound, collapse = "', '"), "'.")
  }
  
  color.vals <- generateColors4Temps(configTable, paletteName = paletteName)
  # extract only those rows which passed the filter criteria
  filterCol <- grep("passed_filter", colnames(dataTable), value = TRUE) # new
  dataTableNew <- filter_(dataTable, filterCol)
  idsFiltered <- unique(dataTableNew[[idVar]])
  
  # loop over all detected proteins and create data.frame for each of them
  plotList <- lapply(idsFiltered, function(prot){
    
    CCR.subset <- dataTableNew[which(dataTableNew[[idVar]] == prot),]
    
    if (nrow(CCR.subset)!=0){
      if (verbose){
        message(paste("Plotting curves for ", prot, sep="")) 
      }
      # process colnames
      colnames(CCR.subset)[grep("transformed", colnames(CCR.subset))] <- 
        sub(paste("(",fcStr, ")",  "(.*)(_transformed)", sep=""), "\\2", 
            colnames(CCR.subset)[grep("transformed", colnames(CCR.subset))])
      
      # create melted data frame
      idVarsForMelting <- c("experiment", "temperature","pEC50", "slope", "R_sq") # new
      m.df <- melt(CCR.subset, id.vars=c(idVar, idVarsForMelting),
                   measure.vars=colnames(CCR.subset)[which(colnames(CCR.subset) %in% 
                                                             unique(extractConc(configTable)))])
      
      # create information column for grouping and annotation of curves
      m.df$info <- sapply(seq(nrow(m.df)), function(row){
        return(paste(paste(m.df$temperature[row], "C", sep=" "),
                     paste("pEC50=", round(m.df$pEC50[row], 2), sep=""),
                     paste("\n", "Hill=", round(m.df$slope[row], 2), sep=""),
                     paste("R2=", round(m.df$R_sq[row], 2), sep=" "), sep=",  "))
      })
      
      # m.df$info <- sapply(seq(nrow(m.df)), function(row){
      #   return(paste(m.df$temperature[row], "C", sep=" "))
      # })
      
      # create df with unique curve fit parameters
      uniq.df <- unique(m.df[which(!is.na(m.df$pEC50)), c("pEC50", "slope", "R_sq", "temperature", 
                                                          "info")])
      
      
      # determine minimal and maximal concentration
      max.conc <- max(log10(10^-6 * as.numeric(as.character(m.df$variable))))
      if (is.infinite(max.conc) | is.na(max.conc)){
        max.conc <- -5
      }
      
      # define x-values for fit
      xvals <- seq(-10, max.conc, 0.01)
      
      # loop over uniq.df and calc sigmoid values for all temperatures
      sigmoid.df <- do.call(rbind, lapply(1:length(uniq.df[,1]), function(i){
        sig.vals <- eval(parse(text=fctSigmoidCCR()), 
                         envir=list(x=xvals,
                                    infl=-as.numeric(as.character(uniq.df$pEC50[i])), 
                                    hill=as.numeric(as.character(uniq.df$slope[i]))))
        temp.df <- data.frame(xvals, sig.vals, temperature=rep(uniq.df$temperature[i], length(xvals)),
                              info=rep(uniq.df$info[i], length(xvals)))
        return(temp.df)
      }))
      
      temp.melt.plot <- ggplot(data=m.df, aes(x=log10(10^-6 * as.numeric(as.character(variable))), 
                                              y=value, colour=as.factor(info))) + 
        geom_point() + 
        ggtitle(paste(m.df[[idVar]][1], sep=" | ")) +
        xlab(paste(compound, "conc. (log M)", sep=" ")) + 
        ylab("Normalized apparent stability") +
        xlim(-10, max.conc) +
        scale_colour_manual(name="Temperature", 
                            values=as.character(color.vals[which(names(color.vals) %in% m.df$temperature)])) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
          axis.line.y = element_line(colour = "black", size=0.5, linetype="solid")) +
        geom_line(data=sigmoid.df, aes(x=xvals, y=sig.vals, colour=as.factor(info)))
      
      
      return(temp.melt.plot)
    } else {
      if (verbose){
        message(paste("None of the curves fulfilled the quality criteria for ", prot, sep=""))
      }
      return(NULL)
    }
  })
  
  names(plotList) <- unique(dataTableNew[[idVar]])
  return(plotList)
}
