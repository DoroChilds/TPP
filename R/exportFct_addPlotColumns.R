exportFct_addPlotColumns <- function(tab, path, idVar, trRef){
  ## Create plots paths to each respective pdf and check for existance.
  
  ## Preparation
  idsProcessed <- tab %>% 
    # Extract id column as character vector:
    extract2(idVar) %>% 
    # Replace special characters by '_':
    gsub("(\\.)", "_", .) 
  
  # Only include reference data plots if a valid reference data set is defined:
  existsRef <- !is.null(trRef) && file.exists(trRef)
  
  ## All plot types share an initial prefix consisting of directory + protein id: 
  pathBase <- file.path(path, "plots", idsProcessed)
  
  ## Create plot paths:
  paths_allCurves    <- paste0(pathBase, "_2D_TPP_all_plots.pdf")
  paths_goodCurves   <- paste0(pathBase, "_2D_TPP_good_plots.pdf")
  paths_singleCurves <- paste0(pathBase, "_2D_TPP_single_plots.pdf")
  paths_splineCurves <- paste0(pathBase, "_2D_TPP_spline_plots.pdf")
  
  ## Check files for existence and remove duplicates:
  paths_allCurves_final     <- removeInvalidPaths(paths_allCurves)
  paths_goodCurves_final    <- removeInvalidPaths(paths_goodCurves)
  paths_singleCurves_final  <- removeInvalidPaths(paths_singleCurves)
  paths_splineCurves_final  <- removeInvalidPaths(paths_splineCurves)
  
  ## Add to output table:
  tab$plot_all_drcurves   <- paths_allCurves_final
  tab$plot_good_drcurve   <- paths_goodCurves_final
  tab$plot_single_drcurve <- paths_singleCurves_final
  tab$plot_spline_fits    <- paths_splineCurves_final
  
  ## Special case: reference data boxplots
  if (existsRef){
    ### Potential problem: in Nils original code, the ids were not preprocessed 
    ###                    by special character removal. Test for consistency 
    ###                    with the plot generating function (still to be written).
    dirRefBoxplots <- file.path(dirname(trRef), "fcBoxplots")
    filesRefBoxplots <- paste0("fcBoxpl_", idsProcessed) 
    paths_refBoxplots <- file.path(dirRefBoxplots, filesRefBoxplots)
    paths_refBoxplots_final <- removeInvalidPaths(paths_refBoxplots)
    tab$plot_tr_reference <- paths_refBoxplots_final
  }
  
  return(tab)
}
