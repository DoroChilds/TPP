#' @title TPP-TR reference object
#'   
#' @description Definition of a TPP-TR reference object
#'  
#' @return A TPP-TR reference object
#' 
#' @examples 
#' trRef <- system.file("example_data/2D_example_data/referenceNormData.RData", package="TPP")
#' tpp2dTRReferenceObject(tppRefDataPath=trRef)
#'
#' @param tppRefData TPP-TR reference object that can be directly passed to the function
#' @param tppRefDataPath character string containing a system path to a RData file containing an
#'   TPP-TR reference object
#' @param fcStr character string indicating which columns contain the fold changes
#' @param qualColName character string indicating which column contain protein 
#'  identification quality measures
#' 
#' @export
tpp2dTRReferenceObject <- function(tppRefData=NULL, tppRefDataPath=NULL, fcStr="norm_rel_fc_protein_", qualColName="qupm"){
  
  # pre-define variables to prevent NOTE by devtools::check()
  variable <- NULL
  Protein_ID <- NULL
  condition <- NULL
  tmp <- NULL
  fc <- NULL
  temps <- NULL
  pEC50 <- NULL
  pseudo_temp <- NULL
  scaledX <- NULL
  x <- NULL
  y <- NULL
  
  thisEnv <- environment()
  
  if(is.null(tppRefDataPath) & is.null(tppRefData)){
    stop("tpptrRefDataPath and tppRefData are both NULL, one of them must be provided.")
  } else if(!is.null(tppRefDataPath) & !is.null(tppRefData)){
    stop("tpptrRefDataPath and tppRefData are both provided, please provide only one.")
  }
  
  if(!is.null(tppRefDataPath)){
    load(tppRefDataPath)
  }
  
  tppCfgTable  <- tppRefData$tppCfgTable
  temperatures <- tppRefData$temperatures
  detailData   <- tppRefData$sumResTable$detail
  summaryData  <- tppRefData$sumResTable$summary
  lblsByTemp   <- tppRefData$lblsByTemp
  lbls <- colnames(temperatures)
  temps <- temperatures[1, lbls]
  
  # create the list used to represent an object for this class
  me <- list(
    
    # define environment where this list is defined so
    thisEnv <- thisEnv,
    
    # define the accessors for the data fields.
    getEnv <- function(){
      return(get("thisEnv",thisEnv))
    },

    getTppCfgTable <- function(){
      return(get("tppCfgTable", thisEnv))
    }, 
    
    getTemperatures <- function(){
      return(get("temperatures", thisEnv))
    }, 
    
    getDetailData <- function(){
      return(get("detailData", thisEnv))
    }, 
    
    getSummaryData <- function(){
      return(get("summaryData", thisEnv))
    } , 
    
    
    createFCBoxPlot <- function(protID){  
      stopifnot(protID %in% detailData$Protein_ID)
      
      protData_detail <- subset(detailData, Protein_ID==protID)
      
      # create dataframe 
      boxPlotData <- do.call(rbind,lapply(as.character(lblsByTemp$lbl), function(l){
        ptrn <- paste(fcStr, l, "_", sep="")
        idx <- grep(ptrn, names(protData_detail))
        lbl <- rep(l, length(idx))
        fc <- unlist(protData_detail[1, idx])
        tmp <- rep(lblsByTemp[which(lblsByTemp$lbl == l), "temp"], length(idx))
        condition <- gsub(ptrn, "",  names(protData_detail)[idx])
        return(data.frame(lbl, tmp, fc, condition))
      }))

      if(sum(!is.na(boxPlotData$fc)) > 0){
        validConds <- unique(boxPlotData$condition[!is.na(boxPlotData$fc)])
        naConds <- unique(boxPlotData$condition[is.na(boxPlotData$fc)])
        maxL <- max(length(validConds), length(naConds))
        conds <- unique(boxPlotData$condition)
        
        for(cnd in conds){
          newRow <- data.frame("valid_FCs" = ifelse(any(!is.na(subset(boxPlotData, condition == cnd)$fc)), "yes", "no"),
                              "NAs" = ifelse(any(is.na(subset(boxPlotData, condition == cnd)$fc)), "yes", "no"),
                              #qualColName = protData_detail[[paste(qualColName, cnd, sep="_")]], 
                              row.names = cnd)
          
          if(exists("myTable")){
            myTable <- rbind(myTable, newRow)
          }else{
            myTable <- newRow
          }
        }
        
        p <- ggplot(boxPlotData, aes(x=tmp, y=fc)) +
              scale_x_discrete(breaks=sort(as.numeric(as.character(lblsByTemp$temp)))) +
              guides(fill="none") +
              geom_boxplot(aes(group=factor(tmp))) +
              ggtitle(protID) + 
              ylab("normalized relative fold change") +
              xlab(paste("temperature [\U00B0", "C]", sep="")) +
              scale_y_continuous(limits=c(0, 1.5)) +
              theme_bw() + 
              annotation_custom(tableGrob(myTable), 
                                #xmin=37,#(max(temps) - min(temps)/2), 
                                #xmax=66.3,##max(temps), 
                                ymin=1.5 - 0.1*maxL, 
                                ymax=1.5)
        
        return(p)
        
      } else {
        stop("No valid fold changes to plot for ", protID,".", sep="")
      }
    },
    
    createMeltPpEC50plot = function(protCCRData=NULL){
      if(!protCCRData$protID %in% detailData$Protein_ID){
        stop(paste("The protein", as.character(protCCRData$protID), "is not present in the 
                   reference data!", sep=" "))
      }
      
      scaleFac <- 2
      wdth <- 0.1
      
      protData_detail <- subset(detailData, Protein_ID==protCCRData$protID)
      
      idx <- grep("meltPoint_", names(protData_detail))
      meltPoint <- unlist(protData_detail[1, idx])
      condition <- gsub("meltPoint_", "", names(protData_detail)[idx])
      densPlotData  <- data.frame(meltPoint, condition)
      
      validConds <- unique(densPlotData$condition[!is.na(densPlotData$meltPoint)])
      naConds <- unique(densPlotData$condition[is.na(densPlotData$meltPoint)])
      plotLims_x_raw <- range(c(densPlotData$meltPoint, 
                               subset(protCCRData$efficacyData, pEC50!=0)[["temp"]]), na.rm=TRUE)
      plotLims_x_adj <- c(floor(plotLims_x_raw[1]-1), ceiling(plotLims_x_raw[2]+1))
      plotLims_y <- c(0, max(9, ceiling(max(protCCRData$efficacyData$pEC50))))
      densPlotData$scaledX <- rep(1, nrow(densPlotData))       
      
      
      protCCRData$efficacyData$pseudo_temp <- protCCRData$efficacyData$temp
      dupesIdx <- duplicated(protCCRData$efficacyData$pseudo_temp)     
      wdthXtra <- wdth*1.3
      if(any(dupesIdx)){
        for(dupe in protCCRData$efficacyData$pseudo_temp[dupesIdx]){
          currDupeIdx <- protCCRData$efficacyData$pseudo_temp == dupe
          dupeCount <- sum(currDupeIdx)
          rng <- (dupeCount-1) * wdthXtra
          newPos <- seq(dupe-(rng/2), dupe+(rng/2), wdthXtra)
          protCCRData$efficacyData$pseudo_temp[currDupeIdx] <- newPos
        }
      }
      
      bp <- ggplot(protCCRData$efficacyData) + 
        geom_bar(aes(x=pseudo_temp, y=pEC50), fill="red", width=wdth, stat="identity", 
                 colour = "red", size = 0.1, alpha=0.4) + 
        scale_x_continuous(limits=plotLims_x_adj, 
                           breaks=sort(unique(c(round(densPlotData$meltPoint, 1), 
                                                protCCRData$efficacyData$temp)))) +
        scale_y_continuous(limits=plotLims_y) + 
        ylab("pEC50") +
        xlab(paste("temperature [\U00B0", "C]", sep="")) +
        annotate("text", label="pEC50s from TPP-CCR scanning", 
                         x=plotLims_x_adj[2], y=plotLims_y[2], hjust=1, vjust=0, size=6, 
                 colour="red" ) +
        annotate("text", label="melting points from TPP-TR reference", 
                         x=plotLims_x_adj[2], y=(plotLims_y[2]-0.5), hjust=1, vjust=0, 
                 size=6, colour="blue" ) +
        ggtitle(as.character(protCCRData$protID))
      
      if(sum(!is.na(densPlotData$meltPoint)) > 1){
        bp <- bp + geom_density(data=subset(densPlotData, !is.na(meltPoint))["meltPoint"], 
                               aes(x=meltPoint),
                               alpha=0.2, 
                               colour="blue", 
                               fill="blue",
                               trim = FALSE,
                               adjust=0.75)
      }
      
      bp <- bp + geom_point(data=subset(densPlotData, !is.na(meltPoint)), 
                           aes(x=meltPoint, y=scaledX), 
                           alpha=0.6, 
                           colour="blue", 
                           fill="blue", 
                           size=7)
      
      return(bp)
    }
    
  )
  
  # define the value of the list within the current environment.
  assign("this",me,envir=thisEnv)
  
  # set the name for the class
  class(me) <- append(class(me),"tpp2dTRReferenceObject")
  return(me)
}