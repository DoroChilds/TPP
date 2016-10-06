writeTRRefOutputTable <- function(tablesList, filePath, fileName, type='xlsx', multiSheet = TRUE){
  
  stopifnot(type %in% c('xlsx', 'txt'))
  
  if(type == 'xlsx'){
    if(multiSheet){
      
      outPath = file.path(filePath, paste(fileName, type, sep='.'))
      
      wb <- createWorkbook()
      
      for(tabName in names(tablesList)){
        addWorksheet(wb, tabName)
        
        linkCols = grep('plot' , colnames(tablesList[[tabName]]), value = TRUE) 
        if(length(linkCols) > 0){
          class(tablesList[[tabName]][[linkCols]]) <- "hyperlink"
        }
        
        writeDataTable(wb, sheet=tabName,  
                       x=tablesList[[tabName]], 
                       startCol=1, startRow=1, 
                       rowNames=FALSE, colNames=TRUE)
      }
      
      tryCatch({
        saveWorkbook(wb, file=outPath, overwrite=TRUE)
        message(paste("OutputTable written to", outPath, " \n"))
      },
      error = function(err){
        message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
        message(err)
      })
      
    } else {
      for(tabName in names(tablesList)){
        outPath = file.path(filePath, paste(fileName, '_', tabName, '.', type, sep=''))
        
        wb <- createWorkbook()
        
        linkCols = grep('plot' , colnames(tablesList[[tabName]]), value = TRUE) 
        if(length(linkCols) > 0){
          class(tablesList[[tabName]][[linkCols]]) <- "hyperlink"
        }
        
        addWorksheet(wb, tabName)
        writeDataTable(wb, sheet=tabName,  
                       x=tablesList[[tabName]], 
                       startCol=1, startRow=1, 
                       rowNames=FALSE, colNames=TRUE)
        
        tryCatch({
          saveWorkbook(wb, file=outPath, overwrite=TRUE)
          message(paste("OutputTable written to", outPath, " \n"))
        },
        error = function(err){
          message("\nCaution! Excel spreasheet could not be produced correctly due to the following error:")
          message(err)
        })
      }
    }
  } else if(type == 'txt'){
    for(tabName in names(tablesList)){
      outPath = file.path(filePath, paste(tabName, '.', type, sep=''))
      tryCatch({
        write.table(x=tablesList[[tabName]], file=outPath, quote= FALSE , 
                    row.names = FALSE, col.names = TRUE, sep='\t')
        message(paste("OutputTable written to", outPath, " \n"))
      },
      error = function(err){
        message("\nCaution! Text file could not be produced correctly due to the following error:")
        message(err)
      })
    }
  }
}


