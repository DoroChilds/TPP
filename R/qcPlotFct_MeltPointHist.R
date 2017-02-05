qcPlotFct_MeltPointHist <- function(resultTab=resultTable, expNames=expNames, 
                                    minR2, expConds){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  resultTable = comparison <- NULL
  
  if (length(expNames) > 1){
    combis <- combn(expNames,2)
    numCompares <- ncol(combis)  
    
    diffDF = data.frame()
    for (j in 1:numCompares){
      n1 <- combis[1,j] 
      n2 <- combis[2,j]
      
      ## Compute melting point differences (ensure direction 'Treatment' -
      ## 'Vehicle', if information about conditions available)
      cond1 <- unname(expConds[expNames==n1])
      cond2 <- unname(expConds[expNames==n2])
      if (identical(cond2, "Treatment") && identical(cond1, "Vehicle")){
        n1_tmp <- n1
        n1 <- n2
        n2 <- n1_tmp
      }
      mpCol_minuend    <- paste("meltPoint", n1, sep="_")
      mpCol_subtrahend <- paste("meltPoint", n2, sep="_")
      
      mp1 <- resultTab[, mpCol_minuend]
      mp2 <- resultTab[, mpCol_subtrahend]
      mpDiff <- mp1 - mp2
      
      #  filter data according to R_sq
      r2_1 <- resultTab[, paste("R_sq", n1, sep="_")]
      r2_2 <- resultTab[, paste("R_sq", n2, sep="_")]
      mpDiff <- mpDiff[which(r2_1 >= 0.8 & r2_2 > 0.8)]
      
      # Compute median and standard deviation of differences
      mdn = round(median(mpDiff, na.rm=TRUE),3)
      q = quantile(mpDiff, probs = c(0.1587, 0.8413), na.rm=TRUE)
      robustSD_left  = round(q[1], 3)
      robustSD_right = round(q[2], 3)
      
      # Create plot annotation:
      titl=paste('median = ', mdn, ' , \n robust SD: left = ', robustSD_left, 
                 ' | right = ', robustSD_right, sep='')
      
      compName <- paste(n1, "-", n2)
      if (length(mpDiff)>0){ # do not produce plot in case no Tm diff exists for 
        # the given combination (for example, because of 
        # poor curve fits).
        
        
        tmpDF = data.frame('diff' = as.numeric(mpDiff), 'comparison'= compName)
        
        p = ggplot(data = tmpDF)
        p = p + stat_bin(aes(x = diff), na.rm = TRUE, colour = 'black', 
                         alpha = 0.3, geom = "bar", binwidth = 0.2)
        xLabels <- seq(-10,10,by = 5)
        p = p + scale_x_continuous(limits=c(-15,15), breaks=xLabels, labels=xLabels)
        p = p + theme(legend.position="bottom")
        xStr1 <- paste('Melting point difference [\U00B0', 'C]', sep='')
        xStr2 <- paste('(',compName,')', sep='')
        xStr3 <- paste('n =', length(mpDiff), ', R_sq >=', minR2)
        p = p + xlab(paste(xStr1, '\n', xStr2, '\n\n', xStr3))
        p = p + ggtitle(titl)
        print(p)
        
        diffDF = rbind(diffDF, tmpDF)
        #       } else{ # to do!
        #         ggplot(data=data.frame(x=c(-1,-1,1,1), y=c(-1,1,-1,1)), aes(x=x,y=y)) + 
        #           geom_blank() + ylab("Minimal slope") + xlab(xLab) +
        #           geom_text(size=10, aes(x=0, y=0, label="No significant Tm shifts shown here\nbecause p-value calculation was not\npossible for any protein.")) 
      }
    }
    
    p = ggplot(data = diffDF)
    p = p + stat_bin(aes(x = diff, colour = comparison, fill = comparison), 
                     na.rm = TRUE, alpha = 0.05, geom = "area", binwidth = 0.2, 
                     position = position_dodge(width = 0))
    p = p + stat_bin(aes(x = diff, colour = comparison), 
                     na.rm = TRUE, geom = "line", binwidth = 0.2, 
                     position = position_dodge(width = 0))
    xLabels <- seq(-10,10,by = 5)
    p = p + scale_x_continuous(limits=c(-15,15), breaks = xLabels, labels = xLabels)
    p = p + theme(legend.position="bottom")
    p = p + guides(colour = guide_legend(nrow = max(1, floor(numCompares/2) ))) 
    p = p + xlab('melting point difference')
    print(p)
  } else {
    message("No melting point difference histograms were produced because they 
            require at least two experiments.")
  }
  }
