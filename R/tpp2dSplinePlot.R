#' @title Fit splines and generate ggplot visualizations
#'   
#' @description Fit splines through TR reference dataset and extrapolates relative 2D-TPP datapoints, 
#'  then compares spline fits of different treatments with non-treatment with an f-test 
#' 
#' @return A list of ggplots which can be accessed via the unique protein ids in the idVar column 
#' 
#' @examples 
#' load(system.file("example_data/2D_example_data/shortData2d.RData", package="TPP"))
#' trRef <- system.file("example_data/2D_example_data/referenceNormData.RData", package="TPP")
#' 
#' @param data_2D result data.frame from a 2D-TPP CCR analysis
#' @param trRef character string of a valid system path to a TPP-TR reference 
#' RData object
#' @param fcStr character string indicating how columns that will contain the actual 
#'   fold change values will be called. The suffix \code{fcStr} will be pasted in front of
#'   the names of the experiments.
#' @param idVar character string indicating name of the columns containing the unique protein 
#'   identifiers in the 2D data set
#' @param refIdVar character string indicating name of the columns containing the unique protein 
#'   identifiers in the reference data set
#' @param refFcStr character string indicating how columns that will contain the fold change
#'   values in the reference data set
#' @param methods vector of character strings that indicate which methods has been used 
#'   for the previous analysis (default: c("doseResponse"), alternative: c("splineFit") or 
#'   c("doseResponse", "splineFit"))
#' @param verbose print description of problems for each protein for which splines fits could 
#'   not be performed
#' 
#' @export
tpp2dSplinePlot <- function(data_2D=NULL, trRef=NULL, fcStr=NULL, idVar=NULL, 
                            refIdVar = "Protein_ID",
                            methods=c("doseResponse", "splineFit"),
                            refFcStr="norm_rel_fc_protein_", verbose=FALSE){  
  
  ## to do: replace redundant parts by the function 'plot_2D_data_on_temperature_range'
  
  if (is.null(idVar)){
    stop("Please specify idVar! A character string matching the column name of unique protein 
         identifiers!")
  }
  
  message("Generating spline plots ...")
  
  # pre-define variables to prevent NOTE by devtools::check()
  passed_filter = tppRefData = variable = condition = value = temperature =  
    uniqueID = tmp = fc = temps = x = y <- NULL
  
  # extract information von TPP-TR reference
  load(trRef)
  detailData <- tppRefData$sumResTable$detail
  lblsByTemp <- tppRefData$lblsByTemp
  
  # loop over all protIDs
  idVec <- data_2D[[idVar]] %>% as.character %>% unique
  resultList <- lapply(idVec, function(protID){
    # subset 2D data
    protData_2D <- data_2D %>% rename_("uniqueID" = idVar) %>% 
      filter(uniqueID == protID)
    #protData_2D <- data_2D[which(as.character(data_2D[[idVar]])==protID),] 
    if ("doseResponse" %in% methods){
      protData_2D <- protData_2D[, c(grep("normalized_to", colnames(data_2D)), 
                                     grep("temperature", colnames(data_2D)))] 
    }else{
      protData_2D <- protData_2D[, c(grep(fcStr, colnames(data_2D)), 
                                     grep("temperature", colnames(data_2D)))] 
    }
    
    protData_detail = detailData %>% rename_("uniqueID" = refIdVar) %>%
      filter(uniqueID == protID)
    #[which(detailData$Protein_ID==protID),] 
    if (length(which(!is.na(protData_detail)))<10){
      if (verbose){
        message(paste("The TR reference dataset does not supply enough data points for", 
                      protID, sep=" ")) 
      }
      return(NULL)
    }else if (!is.null(protData_2D) && nrow(protData_2D)>4){
      # create dataframe 
      labelVec <- as.character(lblsByTemp$lbl)
      plotList <- lapply(labelVec, function(l){
        ptrn = paste(refFcStr, l, "_", sep="")
        idx = grep(ptrn, names(protData_detail))
        lbl = rep(l, length(idx))
        fc = unlist(protData_detail[1, idx])
        tmp = rep(lblsByTemp[which(lblsByTemp$lbl == l), "temp"], length(idx))
        condition = gsub(ptrn, "",  names(protData_detail)[idx])
        dfTmp <- data.frame(lbl, tmp, fc, condition)
        return(dfTmp)
      })
      plotData <- do.call(rbind,plotList)
      
      # remove rows containing NAs
      plotData <- plotData[complete.cases(plotData),]
      
      sfit <- rlm(plotData$fc ~ ns(as.numeric(as.character(tmp)), df=4), data=plotData, maxit=100)
      
      # calc norm vector
      normVec <- predict(sfit, list(tmp=as.numeric(as.character(protData_2D$temperature))))
      
      # normalize data
      rel2d <- do.call(cbind,lapply(protData_2D, function(x){
        return(as.numeric(as.character(x)))
      }))
      
      if ("doseResponse" %in% methods){
        rel2d[,grep("normalized", colnames(rel2d))] <- as.matrix(rel2d[,grep("normalized", colnames(rel2d))]) * normVec
        colnames(rel2d) <- c(paste(sub(paste(fcStr, "(.*)(_normalized_to_lowest_conc)",
                                             sep=""), "\\1", colnames(rel2d)[seq(ncol(rel2d)-1)]), "uM", sep=""), 
                             "temperature")
      }else{
        rel2d[,grep(fcStr, colnames(rel2d))] <- as.matrix(rel2d[,grep(fcStr, colnames(rel2d))]) * normVec
        colnames(rel2d) <- c(paste(sub(paste(fcStr, "(.*)",
                                             sep=""), "\\1", colnames(rel2d)[seq(ncol(rel2d)-1)]), "uM", sep=""), 
                             "temperature")
      }
      
      rel2d.m <- melt(data.frame(rel2d, check.names=FALSE), id.vars=c("temperature"))
      
      if (nrow(rel2d.m[complete.cases(rel2d.m),])>39){
        # define color platte
        nv <- length(levels(rel2d.m$variable))
        colors <- colorRampPalette(c("orange", "red", "midnightblue"))(nv) 
        
        p <- ggplot(plotData, aes(x=as.numeric(as.character(tmp)), y=as.numeric(as.character(fc)))) +
          geom_point(aes(x=as.numeric(as.character(tmp)), as.numeric(as.character(fc)), 
                         shape=condition), color="gray30") + 
          geom_smooth(data=plotData, method=rlm, linetype=2,
                      formula=(y ~ ns(x, df=4)), se=FALSE, color="gray30", size=0.5) +
          geom_point(data=rel2d.m, aes(x=temperature, y=value, color=variable)) + 
          geom_smooth(data=rel2d.m, method=lm, aes(x=temperature, y=value, color=variable),
                      formula=(y ~ ns(x, df=4)), se=FALSE, size=0.5) +
          scale_colour_manual("concentration", values=colors) + 
          #scale_colour_brewer(type="seq", palette="YlGnBu") +
          xlab("Temperature") + 
          ylab("Normalized apparent stability") +
          ggtitle(protID) +
          theme_classic() +
          theme(axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
                axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"))
        
        return(p)
      }else {
        if (verbose){
          message(paste("There were not enough data points for ", as.character(protID), 
                        " to fit useful splines and perform an f-test!", sep="")) 
        }
        return(NA)
      }
    }else {
      if (verbose){
        message(paste("There were not enough data points for", as.character(protID), 
                      "to fit useful splines!", sep=" ")) 
      }
      return(NULL)
    }
  })
  message("Done.")
  names(resultList) <- unique(data_2D[[idVar]])
  return(resultList)
}
