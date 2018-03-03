#Create plots to visualize variation of abs diff probs across chromosomes and individuals

##################################################
#Load data
absdiff <- readRDS("Data/absdiff.rds")
##################################################

for(i in 1:20) {
  matrix <- absdiff[[i]]
  
  #abs diff vectors for lines
  agg <- apply(matrix,1, mean)

  if(i==1) {abs_matrix <- agg}
  else {
   abs_matrix <- rbind(abs_matrix,agg) 
  }
}

################################################
#by chr
abs_chr <- apply(abs_matrix,1,mean)
names(abs_chr) <- seq(1,20)

png("Plots/aggregate_chr.png")
  barplot(abs_chr, ylab = "Mean Abs Diff", xlab = "Chromosome", main="Mean (across individuals) of Mean (across markers) of Abs Diff")
  abline(h=mean(apply(abs_matrix,1,mean)),col="red", ylim=c(0,0.5))
dev.off()

################################################

#by ind (full)
abs_ind <- apply(abs_matrix,2,mean)
names(abs_ind) <- NA

png("Plots/aggregate_ind.png")
  barplot(abs_ind, ylab = "Mean Abs Diff", xlab = "Individual", main="Aggregate Mean Abs Diff by Ind")
  abline(h=mean(apply(abs_matrix,2,mean)),col="red")
dev.off()

#sorted
png("Plots/aggregate_ind_sorted.png")
  barplot(sort(abs_ind), ylab = "Mean Abs Diff", xlab = "Individual", main="Aggregate Mean Abs Diff by Ind")
  abline(h=mean(apply(abs_matrix,2,mean)),col="red")
dev.off()

#png("Plots/aggregate_ind_hist.png")
#  hist(sort(abs_ind), main="Aggregate Mean Abs Diff by Ind", xlab="Mean abs diff")
#dev.off()




