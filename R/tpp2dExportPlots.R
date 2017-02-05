#' @title Export plots for 2D-TPP experiment.
#' @description Exports plots into plots/ directory in the resultPath 
#' 
#' @return None
#' 
#' @details Creates pdf files of the afore created plots by 
#' \code{plot_2D_data_on_temperature_range} or
#' \code{tpp2dCreateDRplots}
#' 
#' @param plotList list of ggplots returned from one of the plotting functions 
#' @param resultPath path for storing results 
#' @param type character string specifying which type of plot is to be exported
tpp2dExportPlots <- function(plotList, resultPath, type="none"){
  # check whether path does exist
  # if (file.exists(resultPath)){ # new
  # define path to store plots
  plotPath <- file.path(resultPath, "plots")
  if (!file.exists(plotPath)){
    # create plots dir if not existing 
    dir.create(file.path(resultPath, "plots"), recursive = TRUE)
  }
  # loop over all plots and generate single pdfs
  invisible(lapply(names(plotList), function(pl){
    if (!is.null(plotList[[pl]]) && !is.na(plotList[[pl]])){
      if (type=="all" | type=="good"){
        savePl <- try(suppressMessages(
          ggsave(plotList[[pl]], 
                 filename=file.path(plotPath, 
                                    paste(gsub("\\|", "_", 
                                               gsub("(\\.)", "_", pl)),"2D_TPP", 
                                          type, "plots.pdf", sep="_"))
          )))
        if(class(savePl) == "try-error"){
          setwd(plotPath)
          suppressMessages(
            ggsave(plotList[[pl]], 
                   filename=paste(gsub("\\|", "_", 
                                       gsub("(\\.)", "_", pl)),"2D_TPP", 
                                  type, "plots.pdf", sep="_")))
        }
        
      }else if (type=="single"){
        if (!is.null(plotList[[pl]])){
          invisible(lapply(names(plotList[[pl]]), function(temp){
            
            savePl <- try(
              suppressMessages(
                ggsave(plotList[[pl]][[temp]], 
                       filename=file.path(
                         plotPath, paste(gsub("\\|", "_", gsub("(\\.)", "_", pl)), 
                                         gsub("(\\.)", "_", temp), 
                                         "2D_TPP", type, "plots.pdf", sep="_")
                       )
                )
              )
            )
            if(class(savePl) == "try-error"){
              setwd(plotPath)
              suppressMessages(
                ggsave(
                  plotList[[pl]], 
                  filename=paste(gsub("\\|", "_", gsub("(\\.)", "_", pl)), 
                                 gsub("(\\.)", "_", temp), 
                                 "2D_TPP", type, "plots.pdf", sep="_")))
            }
            
          }))
        }
      }else if (type=="spline"){
        savePl <- try(
          suppressMessages(
            ggsave(plotList[[pl]], 
                   filename=file.path(
                     plotPath, paste(gsub("\\|", "_", gsub("(\\.)", "_", pl)), 
                                     "2D_TPP", type, 
                                     "plots.pdf", sep="_")))))
        if(class(savePl) == "try-error"){
          setwd(plotPath)
          suppressMessages(
            ggsave(
              plotList[[pl]], 
              filename=paste(gsub("\\|", "_", gsub("(\\.)", "_", pl)), 
                             "2D_TPP", type, 
                             "plots.pdf", sep="_")))
        }
      }else {
        stop("Please specify a valid argument for 'type' ('all', 'good','single' or 'spline)!")
      }
    }
  }))
  # }else{ # new
  #   stop("Please specify a valid argument for 'resultPath'!")
  # }
  return(NULL)
}
