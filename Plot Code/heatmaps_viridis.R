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
    image(t(matrix), col = viridis(256), zlim=c(0,2),
            xlab = "Marker", ylab= "Individual",main=plot_main, scale=NULL,
            labRow = NA, labCol = NA)
  dev.off()
}

#Presentation images
#Chr 5
matrix <- absdiff[[5]]
filename <- paste0("Presentation/figures/heatmap5.png")

png(filename)
  plot_main <- "Chromosome 5"
  onefifthcol <- ncol(matrix)/5; print(onefifthcol)
  image(t(matrix), col = viridis(256), zlim=c(0,2),
        xlab = "Marker", ylab= "Individual",main=plot_main, #scale=NULL,
        axes=FALSE,cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  axis(1, at = seq(0,1,0.1525088), labels=seq(0,ncol(matrix),1000),srt=45,tick=TRUE)
  axis(2, at = seq(0,1,0.2070393), labels=seq(0,nrow(matrix),100),srt=45,tick=TRUE)
dev.off()

#Chr 16
matrix <- absdiff[[16]]
filename <- paste0("Presentation/figures/heatmap16.png")

png(filename)
  plot_main <- "Chromosome 16"
  onefifthcol <- ncol(matrix)/5; print(onefifthcol)
  image(t(matrix), col = viridis(256), zlim=c(0,2),
        xlab = "Marker", ylab= "Individual",main=plot_main, #scale=NULL,
        axes=FALSE,cex.lab=1.1, cex.axis=1.1, cex.main=1.1, cex.sub=1.1)
  axis(1, at = seq(0,1,0.2296211), labels=seq(0,ncol(matrix),1000),srt=45,tick=TRUE)
  axis(2, at = seq(0,1,0.2070393), labels=seq(0,nrow(matrix),100),srt=45,tick=TRUE)
dev.off()

