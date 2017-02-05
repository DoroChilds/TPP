exportFct_trySave <- function(wb, file){
  tryCatch({
    saveWorkbook(wb, file=file, overwrite=TRUE)
    message("File created successfully!\n")
  },
  error = function(err){
    message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
    message(err)
    message(paste("\n\nAlthough the Excel output failed, you can still access the results of the analysis via its return value or as an R object in the results folder.\n"))
  })
  return(TRUE)
}
