merge_cols <- function(data, fun, ...) {
  # Helper function merge columns of matrix into one vector
  #
  # @param data data matrix
  # @param fun function to deal with different values per row
  # @param ... additional arguments to fun
  #
  # @return vector with combined columns
     
     data <- as.matrix(data)
     
     res <- apply(data, 1, function(xx) fun(unique(na.omit(xx)), ...))
     
     res[res == ""] <- NA

     res     
}
