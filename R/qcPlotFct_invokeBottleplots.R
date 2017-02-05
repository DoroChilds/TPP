qcPlotFct_invokeBottleplots <- function(resultTable, compDF){
  
  ## Initialize variables to prevent "no visible binding for global
  ## variable" NOTE by R CMD check:
  x = y <- NULL
  
  alpha = 0.05 # significance level
  
  ## Retrieve experiment names and annotation:
  for(i in 1:nrow(compDF)){
    expNameV <- compDF[i, "refGroup"]
    expNameT <- compDF[i, "testGroup"]
    compName <- compDF[i, "name"]
    nameMpDiff    <- paste("diff_meltP", compName, sep="_")
    nameMinSlope  <- paste("min_slope",compName, sep="_")
    namePVal      <- paste("pVal_adj", compName, sep="_")
    
    xMpDiff <- resultTable[,nameMpDiff]
    xMinSl  <- resultTable[,nameMinSlope]
    xpVals  <- resultTable[,namePVal]
    isHit   <- xpVals <= alpha
    
    if (any(!is.na(isHit))){
      qcPlotFct_Bottleplot(mpDiffs=xMpDiff, minSlopes=xMinSl, isHit=isHit, 
                           strHit=paste("p_adj <=",alpha), 
                           strNoHit=paste("p_adj >",alpha),
                           expName1=expNameT, expName2=expNameV, addHist=TRUE, 
                           yLimVec=c(0, -1.5))
    } else {
      xStr1 <- paste('Melting point difference [\U00B0', 'C]', sep='')
      xStr2 <- paste('(',expNameT,' - ', expNameV,')', sep='')
      xLab <- paste(xStr1, '\n', xStr2, sep='')
      
      pBlank <- ggplot(data=data.frame(x=c(-1,-1,1,1), y=c(-1,1,-1,1)), aes(x=x,y=y)) + 
        geom_blank() + ylab("Minimal slope") + xlab(xLab) +
        geom_text(size=10, aes(x=0, y=0, label="No significant Tm shifts shown\nhere because p-value calculation\nwas not possible for any protein\nin this comparison.")) 
      print(pBlank)
    }
  }
  return(NULL)
}

