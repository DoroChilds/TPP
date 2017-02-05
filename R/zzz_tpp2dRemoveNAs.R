# #' @title Remove NAs in 2D-TPP data
# #'   
# #' @description Removes NAs in the fold change columns of the 2D-TPP data 
# #' 
# #' @return A dataframe without any NA values
# #' 
# #' @param data.list list of data frames that contain the data for the 2D-TPP experiment after 
# #'   computation of fold changes
# #'   
# #' @export

# tpp2dRemoveNAs <- function(data.list){
#   message("Removing NAs...")
#   corrected.data.list <- lapply(names(data.list), function(l.name){
#     corrected.df <- do.call(rbind, lapply(seq(nrow(data.list[[l.name]])), function(row){
#       if (length(which(is.na(data.list[[l.name]][row,])))==0){
#         return(data.list[[l.name]][row,])
#       }
#     }))
#     return(corrected.df)
#   })
#   names(corrected.data.list) <- names(data.list)
#   return(corrected.data.list)
# }
