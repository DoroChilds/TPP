exportFct_colorCodeColumns <- function(wb, sheet, cols, dat){
  if (length(cols)>0){
    # define excel styles for color coding
    darkgreen   <- createStyle(bgFill="darkgreen")
    lightgreen  <- createStyle(bgFill="darkolivegreen3")
    khaki       <- createStyle(bgFill="khaki1")
    lightorange <- createStyle(bgFill="goldenrod1")
    darkorange  <- createStyle(bgFill="darkorange")
    
    validRows <- 2:(nrow(dat)+1)
    conditionalFormatting(wb, sheet = sheet, cols = cols, rows = validRows, 
                          rule = ">=1.5", style = lightgreen)
    conditionalFormatting(wb, sheet = sheet, cols = cols, rows = validRows, 
                          rule = "<=0.5", style = darkorange)
    conditionalFormatting(wb, sheet = sheet, cols = cols, rows = validRows, 
                          rule = "<=0.67", style = lightorange)
    conditionalFormatting(wb, sheet = sheet, cols = cols, rows = validRows, 
                          rule = "<1.5", style = khaki)
    conditionalFormatting(wb, sheet = sheet, cols = cols, rows = validRows, 
                          rule = ">=2", style = darkgreen)
    # info: last rule appears in Excel as first rule
  }
  return(wb)
}
