#Generate Heatmaps for abs diff in geno probs (between r/qtl2 and doqtl) on each chromosome
library(viridis)

##################################################
#Load data
absdiff <- readRDS("Data/absdiff.rds")
##################################################
#Generate heatmap plots for each chromosome
for(i in 1:20) {
  print(i)
  matrix <- absdiff[[i]]
  filename <- paste0("Plots/heatmap_abs_diff/viridis/abs_heatmap_chr_",i,".png")
  png(filename)
  plot_main <- paste0("Chromosome ",i)
  heatmap(t(matrix),Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,2),
          xlab = "Individual", ylab= "Marker",main=plot_main,
          labRow = NA, labCol = NA)
  dev.off()
}

