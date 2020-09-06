#' @title Plot quality control histograms
#'   
#' @description Plots quality control histograms of pEC50 values of reference dataset and
#'   indicates the pEC50 values of the 2D-TPP experiment
#' 
#' @return A pdf with various quality control plots for a specified 2D-TPP data set 
#' 
#' @param configFile data frame or system path to table that specifies important details 
#'  of the 2D-TPP experiment
#' @param resultTable data.frame containing the results of a CCR analysis of 2D-TPP data
#' @param resultPath character string containing a valid system path to which the the qc
#'  plots will be written
#' @param trRef character string with a link to a TPP-TR reference object RData file
#' @param fcStr character string indicating how columns that will contain the actual 
#'  fold change values are called. 
#' @param idVar character string indicating name of the columns containing the unique protein 
#'  identifiers
#' @param qualColName character string indicating which column contain protein 
#'  identification quality measures 
#' 
#' @export
tpp2dPlotQChist <- function(configFile=NULL, resultTable=NULL, resultPath=NULL, trRef=NULL, 
                            fcStr="rel_fc_", idVar="gene_name", qualColName="qupm"){
  # Problem: this function does not return plots, but stores them as a 'side effect'.
  #         In contrast, tpp2dPlotQCpEC50 returns the plots. 
  #         -> Choose one version for consistency.
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  concentration = experiment = intercept = slope = stats = R2 = temperature = 
    stddev = dmso1_vs_dmso2 = marked = experiment = Var1 = Freq = Var2 = 
    median_meltPoint = rev_cumsum <- NULL
  
  message("Creating QC histograms...")  

  # define suffix of normlized data columns
  fcStrNorm <- paste("norm", fcStr, sep="_")
  
  # eval configFile
  configTable <- importCheckConfigTable(infoTable = configFile, type = "2D")
  
  # pre-define global variables
  key <- NULL
  fc <- NULL
  log2fc <- NULL
  type <- NULL
  variable <- NULL
  
  
  # get Experiment IDs from configTable
  expIDs <- unique(configTable$Experiment)
  
  # create directory for histogram qc plots
  dirPath <- file.path(resultPath, "qc_Histograms")
  if (!file.exists(dirPath)){ # new (remove the 'path already exists' warning)
    dir.create(dirPath, recursive = TRUE)
  }
  
  # prepare filename and open pdf for writing
  fileName <- file.path(dirPath, "qc_histograms.pdf")
  pdf(fileName, onefile=TRUE)
  

  temp <- "temperature"
  exp <- "experiment"
  conc <- "concentration"
  type <- "type"
  fc <- "fc"
  log2fc <- "log2fc"
  dmsoRatio <- "dmso1_vs_dmso2"
  extra <- grep(paste0("^", fcStrNorm,".*unmodified"),colnames(resultTable),value=TRUE)
  extra2 <- grep(paste0("^", fcStr),colnames(resultTable),value=TRUE)
  
  # create tidy dataframe
  tmp.df <- resultTable %>%
    select(!!!syms(c(idVar, temp, exp, extra, extra2))) %>%
    gather("key", "fc", !!!syms(c(extra, extra2))) %>%
    mutate(concentration=gsub("(.+)_([0-9,\\.]+)(.*)", "\\2", key),
           type=sub(paste(fcStr, ".*", sep=""), "original", 
                    sub(".*unmodified", "normalized", key)),
           log2fc=log2(fc)) %>%
    select(!!!syms(idVar, temp, exp, conc, type, fc, log2fc)) %>%
    filter(concentration!="0")
  
  # loop over all ms experiments
  lapply(expIDs, function(ms_exp){
    # histogram
    t.plot <- ggplot(tmp.df[which(tmp.df$experiment==ms_exp),], aes(x=log2fc)) + 
      geom_histogram(binwidth=0.1, na.rm = TRUE) + 
      xlim(-2, 2) +
      geom_vline(xintercept=0, color="red", linetype=2) +
      facet_grid(temperature+type~concentration) +
      xlab("log2 FC") + ylab("No. of proteins") +
      ggtitle("original & normalized log2 FC (to resp. 0uM)\nper compound concentration 
              and temperature") + 
      theme_bw() +
      theme(axis.title=element_text(size=20, face="bold"))
    
    grid.arrange(t.plot) 
    
    # density plot
    t.plot <- ggplot(tmp.df %>% filter(experiment==ms_exp), aes(x=log2fc)) + 
      geom_density(adjust=0.1, na.rm = TRUE) + 
      xlim(-2, 2) +
      geom_vline(xintercept=0, color="red", linetype=2) +
      facet_grid(temperature+type~concentration) +
      xlab("log2 FC") + ylab("Density") +
      ggtitle("original & normalized log2 FC (to resp. 0uM)\nper compound concentration 
              and temperature") + 
      theme_bw() +
      theme(axis.title=element_text(size=20, face="bold"))
    
    grid.arrange(t.plot) 
  })
  
   ## test protein FC distributions for normality and create QQ-plot
    tmp.df <- tmp.df[which(tmp.df$type=="normalized"),]
    qq.data <- ddply(tmp.df, c("temperature","concentration"), function(df) {
     # calculate quantiles for FC
     fcq <- quantile(df$log2fc, probs=(1:99)/100, na.rm=TRUE)
     ndq <- qnorm(p=(1:99)/100)
     slope <- (fcq[75]-fcq[25])/(ndq[75]-ndq[25])
     intercept <- fcq[25] - slope*ndq[25]
     
     # calculate R2
     incl <- 1:99 # c(1:25, 75:99)
     ss.res <- sum((fcq[incl] - (ndq[incl]*slope + intercept))^2)
     ss.tot <- sum((fcq[incl] - mean(fcq))^2)
     r2 <- 1 - (ss.res/ss.tot)
     
     # save data in data.frame
     data.frame(slope=slope,
                intercept=intercept,
                R2=round(r2,2),
                stats=paste("R2 =", round(r2,2)))
   })
   
  invisible(lapply(seq(1, length(levels(as.factor(tmp.df$temperature))), by=4), function(i){
    l.temps <- 
      levels(as.factor(tmp.df$temperature))[i:min(i+3,length(levels(as.factor(tmp.df$temperature))))]
    l.qq.tmp <- qq.data[which(qq.data$temperature%in%l.temps),]
    
    # create QQ-plot for FC
    l.plot <- ggplot(tmp.df[which(tmp.df$temperature%in%l.temps),], aes(sample=log2fc)) + 
      facet_grid(temperature~concentration) +
      geom_point(stat="qq", distribution=qnorm, na.rm = TRUE) +
      xlab("Theoret. quantiles from Normal distribution") + ylab("Quantiles of FC") +
      geom_abline(data=l.qq.tmp, aes(intercept=intercept, slope=slope), color="red", linetype=2) +
      xlim(-6,6) + ylim(-6,6) +
      ggtitle("QQ-plot: Normal distribution vs. protein fold changes") + 
      theme_bw() +
      geom_text(data=l.qq.tmp, aes(x=-5, y=5, label=stats), 
                colour="black", inherit.aes=FALSE, parse=FALSE, size=2, hjust=0)

     # save plot
     grid.arrange(l.plot)
  }))     
   
   ## plot R2 of QQ-plot line fit over concetration
   l.plot <- ggplot(qq.data, aes(x=concentration, y=R2, color=as.factor(temperature), 
                                 group=as.factor(temperature))) +
      facet_wrap( ~ temperature, ncol=2) + geom_point() + geom_line() +
      scale_color_discrete("Temperature") +
      xlab("Compound concentration") + ylab("R2 of QQ-plot line fit") +
      geom_hline(yintercept=0.8, colour="red", linetype=2) + 
      #theme(axis.title=element_text(face="bold"), 
    #       axis.text.x=element_text(size=10, face="bold", angle=90), 
    #       axis.text.y=element_text(size=10, face="bold")) +
      theme_bw()+
      ggtitle("R2 of QQ-plot line fit")
   
   # save plot
   grid.arrange(l.plot) 
  
   
   ## plot std. dev. over concentration
   sd.data <- ddply(tmp.df, c("temperature","concentration"), function(df) {
     # calculate std. dev.
     stddev <- sd(df$log2fc, na.rm=TRUE)
     
     # calculate robust sd estimation
     quantiles <- quantile(df$log2fc, c(0.1587, 0.5, 0.8413), na.rm=TRUE)
     mu <- quantiles[2]
     sigma_left <- mu - quantiles[1]
     sigma_right <- quantiles[3] - mu
     
     # save data in data.frame
     data.frame(stddev=round(stddev,2),
                sl=round(sigma_left,2),
                sr=round(sigma_right,2))
   })
   
   # create  plot
  l.plot <- ggplot(sd.data, aes(x=concentration, y=stddev, color=as.factor(temperature), 
                                group=as.factor(temperature))) +
    facet_wrap( ~ temperature, ncol=2) + geom_point() + geom_line() +
    scale_color_discrete("Temperature") +
    xlab("Compound concentration") + ylab("Std. dev. of log2 FC distr.") +
    theme_bw() +
    ggtitle("Std. dev. of log2 FC distribution")
  
  # save plot
  grid.arrange(l.plot) 
  
  if (length(grep(dmsoRatio, colnames(resultTable)))>0){
    tmp.df <- resultTable %>% 
      select(!!!syms(idVar, exp, temp, qualColName, dmsoRatio)) %>%
      filter(!!sym(qualColName)>1) %>%
      mutate(marked=abs(log2(as.numeric(dmso1_vs_dmso2)))>=1,
             annotation=sub("TRUE", "abs(log2(DMSO1/DMSO2))>=1", marked))  %>%
      mutate(annotation=sub("FALSE", "abs(log2(DMSO1/DMSO2))<1", annotation)) %>%
      ddply(.variables=c(idVar, exp, qualColName, dmsoRatio, "marked", "annotation"), 
            function(df) {
              data.frame(temperature=paste(df$temperature, collapse=", "))
            }) %>%
      mutate(label=paste(temperature, "\n(", experiment, ")", sep=""))
    
    tbl <- as.data.frame(table(tmp.df$label, tmp.df$annotation))
    t.plot <- ggplot(tbl) + 
      geom_bar(stat="identity", aes(Var1, Freq, fill=Var2), position="dodge") +
      geom_text(aes(Var1, Freq, label=Freq), position="dodge", vjust=-0.3) +
      scale_fill_manual(values=c("grey","red")) +
      scale_colour_manual(values=c("grey","red")) +
      xlab("MS experiment") + ylab("No. of proteins") +
      ggtitle("No. of proteins with qupm>1") +
      theme_bw() + 
      theme(axis.text.x=element_text(size=12, face="bold", angle=90), 
            axis.title=element_text(size=18, face="bold"))
    
    grid.arrange(t.plot) 
    
    # check if TPP-TR reference is given
    if (!is.null(trRef)) {
      ## DSMO shift plots with reference melting point
      dirName <- rev(strsplit(trRef, "/")[[1]])[2]
      summaryFile <- file.path(sub(basename(trRef), "", trRef), "Summary.txt")
      sumF <- read.table(summaryFile, header=TRUE, sep="\t")
      tmp.df <- merge(tmp.df, sumF[,c("Protein_ID","median_meltPoint")], 
                      by.x=c("representative"), by.y=c("Protein_ID"))
      
      # through out NAs
      tmp.df <- tmp.df[which(complete.cases(tmp.df)),]
      
      # create plot
      t.plot <- ggplot(tmp.df[tmp.df$marked,], aes(x=median_meltPoint)) + 
        geom_histogram(binwidth=0.5) +
        facet_wrap(~label) +
        xlab("Melting temperature") + 
        ylab("No. of proteins") + 
        ggtitle("No. of proteins with qupm>1\n and abs(log2(DMSO1/DMSO2))>=1") + 
        theme_bw() +
        theme(axis.title=element_text(size=20, face="bold"), 
              strip.text.x=element_text(size=14, face="bold"))
      
      # save plot to PDF
      grid.arrange(t.plot)
      
      # create plot
      t.plot <- ggplot(tmp.df[tmp.df$marked,], aes(x=as.factor(experiment), y=median_meltPoint)) + 
        geom_boxplot(aes(fill=experiment)) +
        geom_jitter() + 
        #facet_wrap(~label) +
        ylab("Melting temperature") + 
        xlab("Experiment IDs") + 
        ggtitle("Distribution of proteins with qupm>1\n and abs(log2(DMSO1/DMSO2))>=1") + 
        theme_bw() +
        theme(axis.title=element_text(size=20, face="bold"), 
              strip.text.x=element_text(size=14, face="bold"),
              legend.position=c("none"))
      
      # save plot to PDF
      grid.arrange(t.plot)
    }
  }
  
  # no. of temperatures per protein plot
  tmp.df <- resultTable %>% 
    select(!!!syms(idVar, temp, qualColName)) %>%
    filter(!!sym(qualColName)>1) %>%
    group_by(!!idVar) %>% 
    summarise(count=length(temperature))
  tbl <- as.data.frame(table(tmp.df$count))
  
  # create plot
  t.plot <- ggplot(tbl, aes(x=Var1, y=Freq)) + geom_bar(stat="identity") +
    geom_text(aes(label=Freq), vjust=-0.3, size=6) +
    xlab("No. of temperatures protein was found") + ylab("No. of proteins") +
    ggtitle("No. of proteins with qupm>1") +
    theme_bw() +
    theme(axis.title=element_text(size=20, face="bold"), strip.text.x=element_text(size=12))
  
  # save plot to PDF
  grid.arrange(t.plot)
  # cumulated no. of temperatures per protein plot
  tbl$rev_cumsum <- rev(cumsum(rev(tbl$Freq)))
  
  t.plot <- ggplot(tbl, aes(x=Var1, y=rev_cumsum)) + geom_bar(stat="identity") +
    geom_text(aes(label=rev_cumsum), vjust=-0.3, size=6) +
    xlab("Min. no. of temperatures protein was found") + ylab("No. of proteins") +
    ggtitle("No. of proteins with qupm>1") +
    theme_bw() +
    theme(axis.title=element_text(size=20, face="bold"), strip.text.x=element_text(size=12))
  
  grid.arrange(t.plot)
  
  # close pdf
  dev.off()
  message("Done.")
}
