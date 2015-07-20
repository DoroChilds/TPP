#' @title Fit and analyse dose response curves
#' @description \code{tppccrCurveFit} fits logistic dose response curves to fold
#'   change measurements of a TPP-CCR experiment and returns quality information
#'   about the estimated parameters.
#' @seealso \code{\link{tppccrImport}}, \code{\link{tppDefaultTheme}}
#'   
#' @return A data frame in which the fit results are stored row-wise for each
#'   protein.
#'   
#' @examples
#' data(hdacCCR_smallExample)
#' tppccrData <- tppccrImport(configTable=hdacCCR_config_repl1, 
#'                            data=hdacCCR_data_repl1)
#' tppccrNorm <- tppccrNormalize(data=tppccrData)
#' tppccrTransformed <- tppccrTransform(data=tppccrNorm)
#' tppccrFitResults <- tppccrCurveFit(data=tppccrTransformed)
#' subset(tppccrFitResults, passed_filter)
#' 
#' @param data ExpressionSet containing protein fold changes for curve fitting.
#' @param resultPath location where to store dose-response curve plots and
#'   result table.
#' @param ggplotTheme ggplot theme for dose response curve plots.
#' @param doPlot boolan value indicating whether dose response curves should be
#'   plotted. Deactivating plotting decreases runtime.
#' @param fcCutoff cutoff for highest compound concentration fold change.
#' @param r2Cutoff quality criterion on dose response curve fit.
#' @param slopeBounds bounds on the slope parameter for dose response curve
#'   fitting.
#' @export
#' @details \code{data} is an ExpressionSet object created by
#' \code{\link{tppccrImport}}. It contains the isobaric labels and administered
#' drug concentrations in the \code{phenoData} and user-defined protein
#' properties in the \code{featureData}. Protein IDs are stored in the
#' \code{featureNames}.
#' 
#' The dose response curve plots will be stored in a subfolder with name 
#' \code{DoseResponse_Curves} at the location specified by \code{resultPath}.

tppccrCurveFit <- function(data, resultPath=NULL, ggplotTheme=tppDefaultTheme(), 
                           doPlot=TRUE, fcCutoff=1.5, r2Cutoff=0.8, 
                           slopeBounds=c(1,50)){
  message("Calculating pEC50s ... ")
  if (is.null(resultPath)){
    doPlot <- FALSE
  }
  
  concentrations <- log10(data$concentration * 10^-6)
  concentrations[is.infinite(concentrations)] <- -15
  lbnd <- concentrations[order(concentrations)][2]
  ubnd <- concentrations[order(concentrations)][length(concentrations)]
  label_order  <- data$label
  expName      <- annotation(data)[["name"]]
  fcMat        <- exprs(data)
  cpdEffects   <- featureData(data)$CompoundEffect
  protIDs      <- featureNames(data)
  
  
  # calculate pEC50s  
  if (doPlot){
    ## Check if output directory exists already. If not, create it here.
    plotDir <- "DoseResponse_Curves"
    if (!file.exists(file.path(resultPath, plotDir))){
      dir.create(file.path(resultPath, plotDir), recursive=TRUE)    
    }
    theme_set(ggplotTheme)
    yMax <- 1.5
  } 
  
  ##  Set concentration boundaries: 
  ## [lowest conc. - 0.5 * dilution  step; highest conc. - 0.5 * dilution step]
  conc_bds <- guessDRparamBoundaries(concentrations)  
  
  # Calculate pEC50s & Hill slope for each protein
  t1 <- Sys.time()
  fitResults <- data.frame()
  for(i in 1:nrow(data)) {
    ipi        <- protIDs[i]
    response   <- fcMat[i,]
    cpd_effect <- cpdEffects[i]
    
    pec50 = pec50_final = hill = r2 = file_name = comment <- NA
    
    # remove NA data points
    not_na <- which(!is.na(response))      
    
    # proceed if transformed response values are available
    if(length(not_na)>0) {
      message("Fitting dose response curve for protein: ", ipi)
      
      ## Retrieve dose (independent variable) and response (dependent variable) 
      ## values:
      response <- unname(unlist(response[not_na]))
      dose     <- concentrations[not_na]
      
      # give initial guess on pEC50 & Hill slope
      boundsTmp <- slopeBounds
      if(cpd_effect=="destabilized") {
        boundsTmp <- sort(-1*boundsTmp)
      }
      pec50_init <- guessInitialpEC50(dose, response, conc_bds)
      hill_init  <- guessInitialDRslope(dose, response, boundsTmp, cpd_effect)
      fit <- fitSigmoidCCR(xVec=dose, yVec=response, 
                           hill_init=hill_init, pec50_init=pec50_init,
                           slopeBounds=boundsTmp, concBounds=conc_bds)
      
      if(class(fit) != "try-error") {
        # if fit was successful extract pEC50 and Hill slope and calculate R2
        hill  <- signif(coef(fit)["hill"], 3)
        pec50 <- signif(coef(fit)["infl"], 3)
        r2    <- signif(rSquared(fit, response), 3)
        pec50_final <- -1*pec50
        
        # plot curve
        if (doPlot){
          xlim <- c(concentrations[2], concentrations[length(concentrations)])
          ylim <- c(-0.5,yMax)
          xvals <- seq(xlim[1],xlim[2], 0.01)  
          yvals <- eval(parse(text=fctSigmoidCCR()), 
                        envir=list(x=xvals, infl=pec50, hill=hill))
          ylab <- "normalized apparent stability\n(relative to highest ligand concentration)"
          if(cpd_effect=="destabilized") {
            ylab <- "normalized apparent stability\n(relative to vehicle)"
          }
          p <- ggplot()
          p <- p + scale_y_continuous(limits=ylim) + 
            scale_x_continuous(limits=xlim)
          p <- p + ggtitle(ipi)
          p <- p + xlab("cpd. conc. (log M)") + ylab(ylab)
          p <- p + geom_point(data=data.frame(dose, response), 
                              mapping=aes(x=dose, y=response), na.rm=TRUE, size=4)
          p <- p + geom_line(data=data.frame(xvals, yvals), 
                             mapping=aes(x=xvals, y=yvals), colour="red")
          
          tableDF <- data.frame(condition = expName,
                                pEC50   = signif(pec50_final,2),
                                slope   = signif(hill,2),
                                R2      = signif(r2,2))
          p <- addTableToPlot(plotObj=p, tableDF=tableDF, 
                              meltVar="condition", clrs="red")
          
          
          file_name  <- file.path(plotDir, paste("drCurve_", 
                                                 gsub("([^[:alnum:]])", "_", 
                                                      ipi),".pdf", sep=""))
          
          ## Print plot to PDF:
          pdf(file=file.path(resultPath, file_name), width=7.87, height=9.84, 
              useDingbats=FALSE)
          grid.arrange(p)
          dev.off()      
          ## We used to plot with the ggsave function until a major update in the 
          ## gridExtra package (to version 2.0.0) made ggsave incompatible with 
          ## gridArrange output. We'll later try to switch back to the command:
          # ggsave(filename=file.path(resultPath, file_name), plot=p, width=20, 
          # height=25, units="cm")
        }
      }
    }
    
    ##  Check if protein passed filter criteria
    passed_filter_fc <- !is.na(cpdEffects[i])
    passed_filter_r2 <- !is.na(r2) & r2>=r2Cutoff
    passed_both <- passed_filter_fc & passed_filter_r2
    
    ## Check if slope is outside of concentration range
    flagOutsideConcRange <- -1*pec50_final > ubnd | -1*pec50_final < lbnd
    flagStr <- ifelse(flagOutsideConcRange, yes="Yes", no="No")
    
    ## Save fitting results
    fitResults <- rbind(fitResults, data.frame(Protein_ID=ipi, pEC50=pec50_final,
                                               slope=hill, R_sq=r2, 
                                               passed_filter=passed_both,
                                               pEC50_outside_conc_range=flagStr,
                                               Plot=file_name, 
                                               stringsAsFactors=FALSE))      
  }
  timeDiff <- Sys.time()-t1
  message("Runtime: ", round(timeDiff, 2), " ", units(timeDiff), "\n")
  
  ## ---------------------------------------------------------------------------
  ## prepare output table
  colnames(fcMat) <- paste(colnames(fcMat), "transformed", sep="_")
  
  if (!doPlot){
    fitResults$Plot <- NULL
  }
  
  ## Append fold changes and fit results for the output table:
  dfFCs   <- data.frame(Protein_ID=protIDs, fcMat, stringsAsFactors=FALSE)
  output <- join(dfFCs, fitResults, by="Protein_ID")
  
  ## Append columns imported from input data:
  dfFeatureDat <- pData(featureData(data))
  dfFeatureDat$CompoundEffect <- NULL
  ## Make column names unique:
  colnames(dfFeatureDat) <- paste(colnames(dfFeatureDat), expName, sep="_")
  dfFeatureDat2 <- data.frame(Protein_ID=featureNames(data), dfFeatureDat, 
                              stringsAsFactors=FALSE)
  output <- join(output, dfFeatureDat2, by="Protein_ID")
  
  message("done.")
  return(output)
}
