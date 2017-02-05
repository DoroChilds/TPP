tpptrQCplots <- function(resultTab, expNames, expConditions, compDF, minR2=0.8, 
                         ggplotTheme=tppDefaultTheme()){
  theme_set(ggplotTheme)
  
  ## 1. Illustrate experiment numbers by Venn diagrams:
  message("Creating venn diagrams...")
  qcPlotFct_VennWrapper(resultTab=resultTab, expNames=expNames, 
                        grConditions=expConditions, compDF=compDF, minR2=minR2)
  
  message("done.\n")
  
  ## 2. QC plot to visualize distribution of melting curve parameters:
  if(!is.null(expConditions) && !is.null(compDF)){
    message("Creating QC plots to visualize minimal slope distributions...")
    suppressWarnings(
      qcPlotFct_invokeBottleplots(resultTable=resultTab, compDF=compDF))
    message("done.\n")
  }
  
  ## 3. QC plot of melting point difference histograms between all experiments
  message("Creating QC plots to visualize differences in melting points...")
  qcPlotFct_MeltPointHist(resultTab=resultTab, expNames=expNames, minR2=minR2,
                          expConds=expConditions)
  message("done.\n")
}
