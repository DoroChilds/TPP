# to do...
# plot_splines_to_file <- function(plots, plotPath, highlightIDs, highlightTxt){
#   
#   pageVec <- sort(unique(plots$pageNum))
#   
#   message("Plotting spline fits for ", nrow(plots), 
#           " identifiers to file '", plotPath, "' .")
#   
#   message("Each page will display 20 individual plots (", 
#           length(pageVec)," pages in total).")
#   
#   # Plot to file
#   pdf(file = plotPath, width = 11.811, height = 7.87, useDingbats = FALSE)
#   for (pageTmp in pageVec){
#     
#     message("Plotting page ", pageTmp)
#     
#     plotsTmp <- filter(plots, pageNum == pageTmp) %>%
#       mutate(uniqueID = factor(uniqueID, levels = plots$uniqueID)) %>%
#       group_by(uniqueID) %>%
#       do(plot = {
#         p <- .$plot[[1]]
#         p + theme(legend.position="none")
#       })
#     
#     grid.arrange(grobs = plotsTmp$plot)
#     print(p)
#   }
#   dev.off()
# }


# plot_splines_to_file <- function(data, modelPredH0, modelPredH1, testResPlot, 
#                                  plotPath, highlightIDs, highlightTxt){
#   
#   pageVec <- sort(unique(testResPlot$pageNum))
#   message("Plotting spline fits for ", nrow(testResPlot), " identifiers to file '", plotPath, "' .")
#   message("Each page will display 20 individual plots (", length(pageVec)," pages in total).")
#   
#   # Sort items the same way as they are in 'testResPlot' to keep the same order for plotting 
#   # (i.e. when proteins were sorted by p-value)
#   sortIdx <- testResPlot %>% select(uniqueID) %>% mutate(idx = 1:nrow(.))
#   data2 <- left_join(data, sortIdx, by = "uniqueID") %>% arrange(idx)
#   
#   # Plot to file
#   pdf(file = plotPath, width = 11.811, height = 7.87, useDingbats = FALSE)
#   for (pageTmp in pageVec){
#     message("Plotting page ", pageTmp)
#     testResTmp <- filter(testResPlot, pageNum == pageTmp)
#     idsTmp <- testResTmp$uniqueID %>% as.character
#     datTmp <- filter(data2, uniqueID %in% idsTmp)
#     predH0Tmp <- filter(modelPredH0, uniqueID %in% idsTmp)
#     predH1Tmp <- filter(modelPredH1, uniqueID %in% idsTmp)
#     p <- create_spline_plots(datTmp = datTmp, predH0Tmp = predH0Tmp, 
#                              predH1Tmp = predH1Tmp, testResTmp = testResTmp, 
#                              highlightIDs = highlightIDs, 
#                              highlightTxt = highlightTxt)
#     print(p)
#   }
#   dev.off()
# }
