tpp2dGeneratePlotLinkCols <- function(tab=NULL, path=NULL, idVar=NULL, trRef=NULL){
  # generate plotPath
  plotPath <- file.path(path, "plots")
  # loop over all rows of the data frame and generate character path to respective pdf of all plots
  if (!is.null(trRef) && file.exists(trRef)){
    tab$plot_tr_reference <- sapply(seq(nrow(tab)), function(row){
      fp <- file.path(gsub(basename(trRef), "", trRef), "fcBoxplots" , paste("fcBoxpl",
                                      as.character(tab[[idVar]][row]), sep="_"))
      if (file.exists(fp)){
        return(fp)
      } else{
        return(NA)
      }
    })
    tab$plot_tr_reference[duplicated(tab$plot_tr_reference, imcomparables=NA)] <- NA
  }
  tab$plot_all_drcurves <- sapply(seq(nrow(tab)), function(row){
    fp <- file.path(plotPath, paste(gsub("(\\.)", "_", as.character(tab[[idVar]][row])), 
                                      "2D_TPP_all_plots.pdf", sep="_"))
    if (file.exists(fp)){
      return(fp)
    } else{
      return(NA)
    }
  })
  tab$plot_all_drcurves[duplicated(tab$plot_all_drcurves, imcomparables=NA)] <- NA
  # 
  tab$plot_good_drcurve <- sapply(seq(nrow(tab)), function(row){
    fp <- file.path(plotPath, paste(gsub("(\\.)", "_", as.character(tab[[idVar]][row])), 
                                      "2D_TPP_good_plots.pdf", sep="_"))
    if (file.exists(fp)){
      return(fp)
    } else{
      return(NA)
    }
  })
  tab$plot_good_drcurve[duplicated(tab$plot_good_drcurve, imcomparables=NA)] <- NA
  #
  tab$plot_single_drcurve <- sapply(seq(nrow(tab)), function(row){
    fp <- file.path(plotPath, paste(gsub("(\\.)", "_", as.character(tab[[idVar]][row])), 
                                    gsub("(\\.)", "_", 
                                         as.character(tab[row, grep("temperature", colnames(tab))])),
                                     "2D_TPP_single_plots.pdf", sep="_"))
    if (file.exists(fp)){
      return(fp)
    } else{
      return(NA)
    }
  })
  tab$plot_spline_fits <- sapply(seq(nrow(tab)), function(row){
    fp <- file.path(plotPath, paste(gsub("(\\.)", "_", as.character(tab[[idVar]][row])), 
                                    "2D_TPP_spline_plots.pdf", sep="_"))
    if (file.exists(fp)){
      return(fp)
    } else{
      return(NA)
    }
  })
  tab$plot_spline_fits[duplicated(tab$plot_spline_fits, imcomparables=NA)] <- NA
  
  return(tab)
}