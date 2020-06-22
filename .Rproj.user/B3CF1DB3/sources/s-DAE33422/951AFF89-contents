#' @title plotPredictions
#' @description Plot Time distributions for predicted individuals
#' @details The function creates plots for showing TSD predictions
#' 
#' @param predObj An object returned from either posteriorTime_NaiveMLE or posteriorTime_FullBayesian
#' @param TrueTimes Vector of times (true). Must be same length as the number of samples
#' @param text_main A text used to show in plot
#' @param text_xlab A text used to show in plot
#' @param text_ylab A text used to show in plot
#' @param pdfFileName File name of pdf to print (NULL means no pdf is created)
#' @importFrom grDevices dev.off
#' @importFrom graphics legend lines abline
#' @export

plotPredictions = function(predObj,TrueTimes=NULL,text_main="", text_xlab="time",text_ylab="Probability density function",pdfFileName=NULL) {

  sites = names(predObj$pdf_UNIsite) #obtain name of sites
  IndNames = rownames(predObj$pdf_MULTI) #obtain name of predicting individuals (if any)
  nInds = nrow(predObj$pdf_MULTI) #number of individuals
  time_grid =  as.numeric(colnames(predObj$pdf_MULTI)) #obtain time grid
  maxYaxis = max(predObj$pdf_MULTI, sapply(predObj$pdf_UNIsite,max))
  
  writeToPdf = ifelse(is.null(pdfFileName),FALSE,TRUE)
  if(writeToPdf) pdf(paste0(pdfFileName,".pdf"),width=12,height=7) #print fig to pdf

  #for each individual  
  for(i in 1:nInds) {
    name = paste("Sample",i) #get name of sample (individual)
    if( !is.null(IndNames[i]) ) name = IndNames[i]
    
    plot(0,0,ty="n",xlim=range(time_grid),ylim=c(0,maxYaxis),xlab=text_xlab, ylab=text_ylab,main=text_main)
    for(j in 1:length(sites)) { #for each sites
      pdf = predObj$pdf_UNIsite[[sites[j]]][i,]
      lines(time_grid,pdf,col=j,lty=2)
    }
    pdf = predObj$pdf_MULTI[i,]
    lines(time_grid,pdf,col=1,lty=1,lwd=2)
    if(!is.null(TrueTimes)) {
      if(length(TrueTimes)!=nInds) stop("TrueTimes argument must have length equal the number of samples")
      abline(v=TrueTimes[i],col="gray",lty=2)
    }
    #legend("topleft",legend="Combined",lwd=2,lty=1)
    legend("topright",legend= c("Combined",sites),col=c(1,1:length(sites)),lty=c(1,rep(2,length(sites))),cex=0.8)
  }
  if(writeToPdf) dev.off()
}
  