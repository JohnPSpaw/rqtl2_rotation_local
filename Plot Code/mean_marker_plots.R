##################################################
#Load data
absdiff <- readRDS("Data/absdiff.rds")
##################################################

##################################################
#Plot 
for(i in 1:20) {
  #load data
  matrix <- absdiff[[i]]
  
  #abs diff vectors for lines
  agg <- apply(matrix,2,mean)
  x <- seq(1,length(agg))
  
  if(i==20){i <- "X"}
  main_title <- paste0("Chromosome ",i) 
  pdf_name_marker <- paste0("Plots/agg_marker/agg_marker_",i,".pdf")
  
  pdf(pdf_name_marker)
  
  #Create empty plot
  plot(0, 0, type="n", xlab="Marker", ylab="Sum of Absolute Difference", main=main_title,
       xaxt="s", yaxt="s", ylim=c(0,2), xlim=c(0,length(agg)))
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col="gray90", border=NA)
  
  #gridlines
  abline(h=seq(0,2,0.25), col="white")
  
  #plot abs diff lines
  lines(x,agg, lwd=1.5, col="red")
  
  dev.off()
}
##################################################


