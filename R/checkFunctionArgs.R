checkFunctionArgs <- function(functionCall, expectedArguments){
  ## Check whether all expected arguments are defined when calling a function
  
  myArgs <- names(functionCall)
  
  # Check for missing function arguments (manual checks with 'hasArg' would also work)
  sapply(expectedArguments, function(arg) {
    if (!arg %in% myArgs){
      stop("Error in ", paste(functionCall)[1], 
           ": argument '", arg, "' is missing, with no default", call. = FALSE)
    }
  })
  
}