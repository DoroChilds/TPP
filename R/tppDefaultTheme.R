#' @title Default ggplot theme for melting curve plots.
#' @description Default theme to be passed to the gplots produced by the TPP package.
#' @details Internally, the theme is used as an argument for the function 
#' \code{ggplot2::theme_set} in order specify the appearance of the melting curve plots.
#'
#' The specified plot properties include bold font and increased font size for axis labels and title, as well as a 90 degree angle for y axis labels.
#'
#' @return ggplot theme with default settings for melting plot appearance.
#' @export
#' @examples
#' # Import data:
#' data(hdacTR_smallExample)
#' tpptrData <- tpptrImport(configTable=hdacTR_config, data=hdacTR_data)
#' # Obtain template with default settings:
#' normRequirements <- tpptrDefaultNormReqs()
#' print(normRequirements)
#' # Relax filter on the 10th fold change column for 
#' # normalization set production:
#' normRequirements$fcRequirements[3,3] <- 0.25
#' # Perform normalization:
#' tpptrNorm <- tpptrNormalize(data=tpptrData, normReqs=)

tppDefaultTheme <- function(){
  basesize = 14
  theme_ppt <- theme_bw(base_size=basesize)
  theme_ppt$axis.text.x = ggplot2::element_text(size=basesize+1 , face="bold")
  theme_ppt$axis.text.y = ggplot2::element_text(size=basesize+1, hjust=1, face="bold")
  theme_ppt$axis.title.x = ggplot2::element_text(size=basesize+2)
  theme_ppt$axis.title.y = ggplot2::element_text(size=basesize+2, angle=90)
  theme_ppt$plot.title = ggplot2::element_text(size=basesize+2, hjust=0.5, face="bold")
  return(theme_ppt)
}
