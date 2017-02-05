create_spline_plots <- function(measurements, predictions, colorBy,
                                highlightIDs, highlightTxt){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  uniqueID = testHypothesis = x = y = colorColumn = highlight <- NULL
  
  # Check for missing function arguments
  checkFunctionArgs(match.call(), c("measurements", "predictions", "colorBy",
                                    "highlightIDs", "highlightTxt"))
  
  if (!("uniqueID" %in% colnames(measurements)))
    stop("'measurements' must contain a column called 'uniqueID'")
  
  if (!("uniqueID" %in% colnames(predictions)))
    stop("'predictions' must contain a column called 'uniqueID'")
  
  ## Sort ID levels so that plots are displayed in the same order as 
  ## they appear in the 'uniqueID' column after separating them by 'facet_wrap':
  measurements <- measurements %>% 
    mutate(uniqueID = factor(uniqueID, levels = unique(uniqueID)))
  
  fctrsH0 <- colorBy %>% 
    filter(testHypothesis == "null") %>% 
    extract2("factors") %>% unique %>% as.character
  
  fctrsH1 <- colorBy %>% 
    filter(testHypothesis == "alternative") %>% 
    extract2("factors") %>% unique %>% as.character
  
  if (length(fctrsH1) == 0){
    fctrsH1 = "condition"
  }
  
  ## Create a column that will be used for coloring by the plot asthetics.
  ## It contains all possible combinations of the factors that distinguish null
  ## and alternative models
    measurements$colorColumn <- subset(measurements, select = fctrsH1) %>%
      apply(., 1, function(x) {
        out <- ifelse(all(is.na(x)), NA, paste(unique(x), collapse = " / "))
      }) %>%
      factor()

  ## Ceck if the predicted values are suitable for plotting. 
  ## Criteria: presence of the column which will be used for coloring 
  ## (named 'colorColumn') and presence of non-NA values.
  ## The color columns are automatically generated from the information in the 
  ## fitted model and therefore missing in case the model fit failed:
  validPredictions <- any(!is.na(predictions$y))
  
  if (validPredictions){
    
    predictions <- predictions %>%
      #filter(!is.na(y)) %>%
      mutate(uniqueID = factor(uniqueID, levels = unique(uniqueID)))
    
    predictionHasAllFactors <- all(fctrsH1 %in% colnames(predictions))
    
    if (!predictionHasAllFactors){
      
      isNullModel <- predictions$testHypothesis == "null"
      
      ## Assign missing column for coloring. Usually, this is done 
      ## by the prediction function which automatically retrieves the factors
      ## for the alternative model fit from the model. In the default use case,
      ## these factors comprise the variable given by fctrsH1. However, they will
      ## not be detected/assigned if only null models were present.
      for (fctr in fctrsH1){
        if (!fctr %in% colnames(predictions)){
          predictions[, fctr] <- ifelse(isNullModel, "null model", "smoothing spline")
        }
      }
    }
    
    ## Create the same color column as in measurements:
    predictions$colorColumn <- subset(predictions, select = fctrsH1) %>%
      apply(., 1, function(x) {
        out <- ifelse(all(is.na(x)), NA, paste(unique(x), collapse = " / "))
        }) %>%
      factor()
    
  }
  
  ## Generate plot scaffold:
  placeholder <- data.frame(x = numeric(), 
                            y = numeric(), 
                            uniqueID = factor(),
                            colorColumn = factor())

  p <- ggplot(data = placeholder, 
              aes(y = y, x = x,  color = colorColumn)) +
    coord_cartesian(ylim = c(-0.1, 2)) +
    ylab("Fraction non-denatured") + 
    xlab("Temperature [\U00B0 C]") +
    scale_shape_discrete(name = "")  +
    theme_bw()
  
  ## Sort colorColumn levels so that a potential 'null model' entry appears
  ## last and gets assigned the black color. This entry was automatically
  ## created by the function 'predict_and_plot_spline_model':
  if(validPredictions){
    
    colorLevels <- c(as.character(measurements$colorColumn), 
                     as.character(predictions$colorColumn))
    
  } else {
    
    colorLevels <- as.character(measurements$colorColumn)
    
  }
  
  colorLevels <- colorLevels %>% unique %>% setdiff(NA)
  
  newColors <- plotColors_splineFits(expConditions = colorLevels)
  newLevels <- names(newColors)
  
  p <- p + scale_color_manual(name = "", values = newColors)
  
  measurements$colorColumn <- factor(measurements$colorColumn, levels = newLevels)
  
  if (validPredictions){
    predictions$colorColumn <- factor(predictions$colorColumn, levels = newLevels)
  }
  
  
  ## Add layers for measurements and predicted curves:
  p <- p + geom_point(data = measurements, 
                      aes(shape = factor(replicate)), na.rm = TRUE)
  
  if (validPredictions){
    p <- p + geom_line(data = predictions, na.rm = TRUE)
  }
  
  ## Check if only a single protein is plotted. If not, split into facets:
  ids <- measurements$uniqueID %>% as.character %>% unique
  
  if (length(ids) >1) {
    p <- p + facet_wrap(~ uniqueID, nrow = 4, ncol = 5) +
      theme(strip.text = element_text(size = 7))
  } else{
    p <- p + ggtitle(ids)
  }
  
  # try(p <- p + geom_label(data = testResTmp, vjust = "top", hjust = "right", 
  #                         label.size = 0.5, inherit.aes = FALSE,
  #                         aes(label = textStr, x = Inf, y = Inf), alpha = 0.1), 
  #     silent = TRUE)
  
  ## Highlight background, if desired (to do: shift this to an extra annotation function, also add the text labels there)
  backgroundTable <- measurements %>% 
    ungroup %>%
    distinct(uniqueID) %>% 
    mutate(highlight = uniqueID %in% highlightIDs)
  
  if (any(backgroundTable$highlight)){
    p <- p + geom_rect(data = backgroundTable, aes(fill = factor(highlight)), 
                       xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                       alpha = 0.1, inherit.aes = FALSE) +
      scale_fill_manual(highlightTxt, values = c("FALSE" = "white", "TRUE" = "green"))
  }
  
  return(p)
}