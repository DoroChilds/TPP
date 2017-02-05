plotIndividual <- function(data, fittedModels, plotAnnotation, 
                           plotDir, filePrefix, returnPlots, nCores){
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("data", "fittedModels", "plotAnnotation", 
                                    "plotDir", "filePrefix", "returnPlots", 
                                    "nCores"))
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  idTmp = uniqueID = textString <- NULL
  
  ## If export to files is desired, define paths first
  doPlot <- !is.null(plotDir)
  
  if (doPlot){
    
    newDir <- file.path(plotDir, "Individual")
    
    if (!file.exists(newDir)){
      dir.create(newDir, recursive=TRUE)
    }
    
    paths <- createValidPlotPaths(type = filePrefix, 
                                  protIDs = unique(data$uniqueID),
                                  plotDir = newDir) %>%
      data.table(key = "uniqueID")
    
  }
  
  
  ## Reformat to data.table in order to speed up subsetting
  data2 <- data.table(data, key = "uniqueID")
  
  fittedModels <- data.table(fittedModels, key = "uniqueID")
  
  if (!is.null(plotAnnotation)){
    plotAnnotation <- data.table(plotAnnotation, key = "uniqueID")
  }
  
  ## Generate plots per protein
  nCores <- checkCPUs(nCores)
  doParallel::registerDoParallel(cores = nCores)
  
  allIDs <- unique(data$uniqueID)
  
  results <- foreach(idTmp = allIDs, .combine=rbind, .inorder=FALSE, 
                     .verbose=FALSE) %dopar% {
                       
                       datTmp <- data2[.(idTmp)]
                       
                       fitsTmp <- fittedModels[.(idTmp)]
                       
                       plotTmp <- predict_and_plot_spline_models(dat = datTmp, fits = fitsTmp)
                       
                       if (!is.null(plotAnnotation)){
                         
                         ## Annotate each plot by p-values and degrees of freedom
                         
                         labelData <- plotAnnotation[.(idTmp)]
                         
                         plotTmp <- plotTmp + 
                           geom_label(data = labelData, 
                                      na.rm = TRUE,
                                      vjust = "top", 
                                      hjust = "right",
                                      label.size = 0.5, 
                                      inherit.aes = FALSE,
                                      aes(label = textString, 
                                          x = Inf, y = Inf), 
                                      alpha = 0.1)
                       }
                       
                       out <- data.frame(uniqueID = idTmp, stringsAsFactors = FALSE)
                       
                       if (returnPlots){
                         out <- out %>% 
                           group_by(uniqueID) %>% 
                           do(plot = plotTmp)
                       }
                       
                       if (doPlot){
                         
                         ## Print plot to PDF
                         
                         fTmp <- paths[.(idTmp)] %>% extract2("path")
                         
                         pdf(file = fTmp, width = 7.87, height = 5.90551, useDingbats = FALSE)
                         print(plotTmp)
                         dev.off()
                         
                         out$path = fTmp
                         
                       }
                       
                       
                       out <- data.table(out)
                       
                       return(out)
                     }
  
  return(results)
}