#' @title Plot single 2D-TPP CCR curves
#' @description Generates a list of plots of all proteins for a certain 2D-TPP experiment
#' 
#' @return A list of all single dose-reponse plots that could be fitted and fullfilled the
#'  requested quality criteria 
#'  
#' @examples
#' data("panobinostat_2DTPP_smallExample")
#' cfg <- panobinostat_2DTPP_config
#' datRaw <- panobinostat_2DTPP_data
#' data2d <- tpp2dImport(cfg, datRaw, fcStr = NULL)
#' fcData2d <- tpp2dComputeFoldChanges(cfg, data2d, intensityStr="sumionarea_protein_")
#' normData2d <- tpp2dNormalize(cfg, fcData2d)
#' config_ccr <- tpp2dCreateCCRConfigFile(cfg)
#' ccr2dResults <- tpp2dCurveFit(config_ccr, normData2d, idVar = "unique_ID")
#' singleCurves <- tpp2dPlotCCRSingleCurves(cfg, ccr2dResults, 
#'                                      idVar = "representative",
#'                                      fcStr = "norm_rel_fc_protein_")
#' singleCurves[["IPI00289601.10"]][[1]] 
#' 
#' @param configTable data frame that specifies important details of the 2D-TPP experiment
#' @param data ouput table returned by the \code{tpp2dCurveFit} function
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param verbose boolean variable stating whether a print description of problems/success for 
#'  plotting of each protein should be printed
#' @export
tpp2dPlotCCRSingleCurves <- function(configTable=NULL, data=NULL, idVar=NULL,
                                     fcStr="rel_fc_protein_", verbose=FALSE){

  # pre-define global variables
  variable <- NULL
  value <- NULL
  sig.vals <- NULL
  
  compound <- configTable$Compound %>% unique
  ## remove the compound name from the result table columns. To avoid confusion, 
  ## throw an error if more than one compound names were in use # new
  if (length(compound) == 1){
    oldCols <- colnames(data)
    patternToRemove <- paste0("_", gsub("[^[:alnum:]]", "_", compound))
    newCols <- gsub(patternToRemove, "", oldCols)
    colnames(data) <- newCols
  } else {
    stop("Can currently only handle config tables with one compound name.
         Compound names detected in your config table: '", 
         paste(compound, collapse = "', '"), "'.")
  }
  
  # loop over all detected proteins and create data.frame for each of them
  singlePlotList <- lapply(unique(data[[idVar]]), function(prot){
    
    CCRsamp <- data[which(data[[idVar]] == prot),]
    
    # extract only those rows which passed the filter criteria
    CCR.subset <- CCRsamp[which(unlist(CCRsamp["passed_filter"])),]
    
    if (nrow(CCR.subset)!=0){
      if (verbose){
        message(paste("Plotting single good curves for ", prot, sep="")) 
      }
      # process colnames
      colnames(CCR.subset)[grep("transformed", colnames(CCR.subset))] <- 
        sub(paste("(",fcStr, ")",  "(.*)(_transformed)", sep=""), "\\2", 
            colnames(CCR.subset)[grep("transformed", colnames(CCR.subset))])
      
      # create melted data frame
      m.df <- melt(CCR.subset, id.vars=c(idVar, "experiment", "temperature","pEC50",
                                         "slope", "R_sq"),
                   measure.vars=colnames(CCR.subset)[which(colnames(CCR.subset) %in% 
                                                             unique(extractConc(configTable)))])
      
      # create df with unique curve fit parameters
      uniq.df <- unique(m.df[which(!is.na(m.df$pEC50)), c("pEC50", "slope", "R_sq", "temperature")])
      
      
      # determine minimal and maximal concentration
      max.conc <- max(log10(10^-6 * as.numeric(as.character(m.df$variable))))
      if (is.infinite(max.conc) | is.na(max.conc)){
        max.conc <- -5
      }
      # define x-values for fit
      xvals <- seq(-10, max.conc, 0.01)
      
      # loop over uniq.df and calc sigmoid values for all temperatures
      sigmoid.df <- do.call(rbind, lapply(seq(nrow(uniq.df)), function(i){
        sig.vals <- eval(parse(text=fctSigmoidCCR()), 
                         envir=list(x=xvals,
                                    infl=-as.numeric(as.character(uniq.df$pEC50[i])), 
                                    hill=as.numeric(as.character(uniq.df$slope[i]))))
        temp.df <- data.frame(xvals, sig.vals, temperature=rep(uniq.df$temperature[i], length(xvals)))
        return(temp.df)
      }))
      
      single.plots <- lapply(unique(m.df$temperature), function(temp){
        sub.m.df <- m.df[which(m.df$temperature == temp),]
        sub.sig.df <- sigmoid.df[which(sigmoid.df$temperature == temp),]
        
        temp.melt.plot <- ggplot(data=sub.m.df, aes(x=log10(10^-6 * as.numeric(as.character(variable))), 
                                                y=value)) + 
          geom_point() + 
          ggtitle(paste(paste(m.df[[idVar]][1], sep=" | "), 
                        paste("T=", temp, "C", sep=" "), 
                        paste("pEC50=", round(sub.m.df$pEC50[1], 2), sep=""), 
                        paste("\n", "Hill=", round(sub.m.df$slope[1], 2), sep=""), 
                        paste("R2=", round(sub.m.df$R_sq[1], 2), sep=" "), sep=",  ")) +
          xlab(paste(compound, "conc. (log M)", sep=" ")) + 
          ylab("Normalized apparent stability") +
          xlim(-10, max.conc) +
          scale_fill_discrete(name="Temperature") +
          theme_classic() +
          theme(
            axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
            axis.line.y = element_line(colour = "black", size=0.5, linetype="solid")) +
          geom_line(data=sub.sig.df, aes(x=xvals, y=sig.vals))
        
        return(temp.melt.plot)
      })
      names(single.plots) <- as.character(unique(m.df$temperature))
      return(single.plots)
      
    } else {
      if (verbose){
        message(paste("None of the curves fullfilled the quality criteria for ", prot, sep="")) 
      }
      return(NULL)
    }
  })
  
  names(singlePlotList) <- unique(data[[idVar]])
  return(singlePlotList)
}
