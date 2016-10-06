#' @title Create data frame list for 2D-TPP experiment
#' @description Creates a 2D-TPP data frame-list featuring a dataframe for each 
#'   temperature analyzed in the experiment.
#'  
#' @return A list of data frames in which the experimental data are stored row-wise for each
#'   protein.
#'   
#' @details Invokes the following steps: \enumerate{ \item  either reads in temperature specific 
#'   experimental data and creates data frame or extracts temperature specific data from pre-
#'   existing data frame list}
#'   
#' @param configTable data frame of experimental conditons. See Section \code{details} for 
#'   instructions how to create this object. 
#' @param data can be either a list of data frames featuring a data frame with data
#'   for each experiment or can be NULL when filepaths for the respective experiments are indicated in 
#'   the configTable
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param intensityStr character string indicating which columns contain the sumionarea values. 
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param addCol additional column names that specify columns in the input data that are 
#'   to be attached to the data frame throughout the analysis 
#'   
#' @export
tpp2dCreateDataFrameList <- function(configTable=NULL, data=NULL, 
                                     idVar="representative", 
                                     fcStr=NULL, addCol=NULL,
                                     intensityStr="sumionarea_protein_", 
                                     qualColName=c("qupm","qusm")){
  # pre-define output list
  #df.list <- list()
  exp.ind <- which(colnames(configTable) %in% "Experiment")
  
  # check whether data is provided or must be read in
  if (is.null(data)){
    # loop over all rows of the config table and read in data corresponding to experimental conditons
    out.df <- lapply(seq(nrow(configTable)), function(exp){
      # get path string
      path.str <- as.character(configTable$Path[exp])
      # get temperature value
      temp.val <- configTable$Temperature[exp]
      # get experiment id
      exp.val <- configTable$Experiment[exp]
      # get columns which correspond to the extracted temperature
      val.cols <- which(configTable[exp,]!="-")[which(!(which(configTable[exp,]!="-") %in% 
                   (which(colnames(configTable[exp,]) %in% 
                     c("Compound", "Experiment", "Temperature", "RefCol", "Path")))))]
      # concatenate fold change value with pre-str
      if (!is.null(fcStr)){
        relevant.cols <- c(idVar, qualColName, addCol,
                           paste(intensityStr, colnames(configTable[,val.cols]), sep=""),
                           paste(fcStr, colnames(configTable[,val.cols]), sep=""))
      }else {
        relevant.cols <- c(idVar, qualColName, addCol,
                           paste(intensityStr, colnames(configTable[,val.cols]), sep=""))
      }
      # read in data frame of respective experiment
      data.f <- importFct_2Ddataframe(filePath=path.str, rowNumber=exp)
      
      # create unique identifier
      u.id <- paste(as.character(configTable[exp, exp.ind]), paste(as.character(temp.val),
                                                                   as.character(data.f[[idVar]]), sep="_"), sep="_")
      
      # throw error if any of the relevant column names can not be found in the column names of 
      # the read-in data frame data.f
      if (!all(relevant.cols %in% colnames(data.f))){
        stop("Some of the columns the import function tried to find in your data do not exist! 
             Please check the suffices and the additional column names you have specified!")
      }
      
      # subset dataframe 
      sub.df <- data.frame(data.f[relevant.cols],
                           temperature=rep(temp.val, nrow(data.f)), 
                           experiment=rep(exp.val, nrow(data.f)),
                           unique_ID=u.id)   
      return(sub.df)
    })
  }else if (!is.data.frame(data[[1]])){
    stop("Please hand over the experimental data as a list of data frames or specifiy respective 
         file paths in the config table!")
  } else {
    # loop over all rows of the config table
    out.df <- lapply(seq(nrow(configTable)), function(exp){
      # get temperature value
      temp.val <- configTable$Temperature[exp]
      # get experiment id
      exp.val <- configTable$Experiment[exp]
      # get columns which correspond to the extracted temperature
      val.cols <- which(configTable[exp,]!="-")[which(!(which(configTable[exp,]!="-") %in% 
                   (which(colnames(configTable[exp,]) %in% 
                    c("Compound", "Experiment", "Temperature", "RefCol", "Path")))))]
      # concatenate fold change value with pre-str
      if (!is.null(fcStr)){
        relevant.cols <- c(idVar, qualColName, addCol,
                           paste(intensityStr, colnames(configTable[,val.cols]), sep=""),
                           paste(fcStr, colnames(configTable[,val.cols]), sep=""))
      }else {
        relevant.cols <- c(idVar, qualColName, addCol,
                           paste(intensityStr, colnames(configTable[,val.cols]), sep=""))
      }
      # create unique identifier
      u.id <- paste(as.character(configTable[exp, exp.ind]), paste(as.character(temp.val),
              as.character(data[[configTable[exp, exp.ind]]][[idVar]]), sep="_"), sep="_")
      # subset dataframe 
      sub.df <- data.frame(data[[configTable[exp, exp.ind] ]][relevant.cols],
                           temperature=rep(temp.val, length(data[[configTable[exp,exp.ind]]][,1])), 
                           experiment=rep(exp.val, length(data[[configTable[exp,exp.ind]]][,1])),
                           unique_ID=u.id)    
      return(sub.df)
    })
  }
  # create data.list names
  dt.names <- sapply(seq(nrow(configTable)), function(exp){
    return(paste(configTable$Experiment[exp], configTable$Temperature[exp], sep="_")) 
  })
  names(out.df) <- dt.names
  return(out.df)  
}