qcPlotFct_VennWrapper <- function(resultTab, expNames, grConditions, compDF, 
                                  minR2){  
  ## Invoke Venn diagram computation to compare protein numbers per experiment,
  ## as well as their overlap.
  
  # supress useless VD logfiles
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  ## --------------------------------------------------------------------------
  ## I) Create a Venn diagram of protein numbers in all experiments
  ## --------------------------------------------------------------------------
  if(length(expNames) <= 5){
    message("Creating Venn diagram to summarize all experiments.")
    mainTxt = paste(nrow(resultTab), 'proteins identified overall')
    plotObj <- qcPlotFct_CreateVenn(dataTab=resultTab, grNames=expNames,
                                    mainTxt=mainTxt, method="individual", 
                                    minR2=minR2)
    grid.newpage()
    grid.draw(plotObj)
  } else {
    message("Venn diagrams can only be created for up to 5 experiments.")
  }
  ## --------------------------------------------------------------------------
  ## II) Creates separate Venn diagrams to compare the experiments in each 
  ## user-defined comparison pair (if applicable)
  ## --------------------------------------------------------------------------
  if (!is.null(compDF)){
    if(nrow(compDF) <= 5){
      message("Creating Venn diagram to summarize all comparisons.")
      plotObj <- qcPlotFct_CreateVenn(dataTab=resultTab, compDF=compDF, 
                                      mainTxt="", method="comparison")  
      grid.newpage()
      grid.draw(plotObj)
    } else {
      message("Venn diagrams can only be created for up to 5 comparisons.")
    }
  }
}

