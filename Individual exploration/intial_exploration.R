#Exploring chr=16, ind=10, markers=(577,870)
#Load data
absdiff <- readRDS("Data/absdiff.rds"); attieDO <- readRDS("Data/base/attieDO_cross.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]; ind_index <- seq(1,dim(absdiff[[1]])[1])

ind <-10; chr <- 16; markers <- c(577, 870)

#convert index of individual to individual name
#ind is index with respect to absdiff
ind_name <- ind_names[ind] #name of ind... such as "DO-101"... primary key between absdiff and attieDO
ind_attieDO <- which(dimnames(attieDO$cross_info)[[1]] == ind_name) #index in of ind in attieDO

#raw genotypes (matrix)
g <-attieDO$geno[[chr]]; fg <-attieDO$founder_geno[[chr]] #individuals x markers ..... 1 = AA, 3 = BB, 2 = AB

#markers
markers <- colnames(g)[markers[1]:markers[2]]

#genotypes 
gg <- g[ind_attieDO, markers]; gg[gg==0] <- NA
fgg <- fg[,markers]; fgg[fgg==0] <- NA

#Proportion of markers matching each founder
proportion_match <- rep(NA,8)
for(i in 1:8) {
  proportion_match[i] <- mean(gg == fgg[i,], na.rm=TRUE) 
}

#RQTL2 Predictions





