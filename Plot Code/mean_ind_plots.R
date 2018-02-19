##################################################
#Load data
absdiff <- readRDS("Data/absdiff.rds")
##################################################
#Plot 
for(i in 1:20) {
  matrix <- absdiff[[i]]
  
  #abs diff vectors for lines
  agg <- apply(matrix,1, mean)
  names(agg) <- NA
  x <- seq(1,dim(matrix)[1])
  
  main_title <- paste0("Chromosome ",i) 
  png_name_ind <- paste0("Plots/agg_ind/agg_ind_",i,".png")
  
  png(png_name_ind)
    hist(agg, main = main_title, breaks=30, 
         xlab="Mean across markers", ylab = "Frequency (individuals)", ylim=c(0,80), xlim=c(0,2))
  dev.off()
}

