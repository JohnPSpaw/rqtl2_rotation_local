##Collection of functions for comparing proportions of matching genotypes
##between individual and each founder and a particular individual


#Load data
absdiff <- readRDS("Data/absdiff.rds")
attieDO <- readRDS("Data/base/attieDO_cross.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])

#generate column names
lcount <- 1
geno_names <- rep(NA,36)
for(i in 1:8) {
  for(j in 1:8) {
    if(i <= j) {
      geno_names[lcount] <- paste0(letters[i],letters[j])
      lcount <- lcount + 1
    }
  }
}

#computes expected genotype at a marker given two specified founders
expected_geno <- function(founder1, founder2)
  {
    result <- founder1
    result[founder1==0 | is.na(founder1) | founder2==0 | is.na(founder2)] <- NA
    result[founder1==1 & !is.na(founder1) & founder2==1 & !is.na(founder2)] <- 1
    result[founder1==3 & !is.na(founder1) & founder2==3 & !is.na(founder2)] <- 3
    result[founder1==1 & !is.na(founder1) & founder2==3 & !is.na(founder2)] <- 2
    result[founder1==3 & !is.na(founder1) & founder2==1 & !is.na(founder2)] <- 2
    
    result
  }

#Output proportion of matchs between individual and each founder combination along a marker region
matching_genos <- function(ind, chr, markers) {
  #convert index of individual to individual name
  #ind is index with respect to absdiff
  ind_name <- ind_names[ind] #name of ind... such as "DO-101"... primary key between absdiff and attieDO
  ind_attieDO <- which(dimnames(attieDO$cross_info)[[1]] == ind_name) #index in of ind in attieDO
  
  #individual raw genotypes
  g <-attieDO$geno[[chr]] #matrix individuals x markers ....... 1 = AA, 3 = BB, 2 = AB
  #founder raw genotypes
  fg <-attieDO$founder_geno[[chr]] #matrix founders x markers ......... 1 = AA, 3 = BB, 2 = AB
  
  #markers
  markers <- colnames(g)[markers[1]:markers[2]]
  #genotypes
  gg <- g[ind_attieDO, markers]
  fgg <- fg[,markers]
  # problem: missing values are coded as 0
  gg[gg==0] <- NA
  fgg[fgg==0] <- NA
  
  count <- 1 
  d <- rep(NA,36)
  for(i in 1:8) {
    #print(paste0("i=",i))
    for(j in 1:8) {
      #print(paste0("j=",j))
      if(i <= j) {
        #print(paste0(i,j))
        d[count] <- mean(gg == expected_geno(fgg[i,], fgg[j,]), na.rm=TRUE)
        count <- count + 1
      }
    }
  }
  
  d
}

explore_genos <- function(chr,ind,markers) {
  #Generate matrix (ind x genotype) of proportion of matches
  print("Generating matches")
  for(i in explore_ind) {
    if(i==explore_ind[1]){
      matching <- round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2)
    }  
    else {
      matching <- rbind(matching, round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2))
    }
  }
  colnames(matching) <- geno_names
  #print(head(matching))
  print(matching[1:20,])
  
  print("Plotting")
  directory <- "Plots/explore_matches/"
  fname <- paste0(directory,"explore_chr",chr,"_markers",markers[1],"-",markers[2],"_inds",ind[1],"-",ind[length(ind)],".png")
  png(fname)
  heatmap(t(matching),Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,1),
          xlab = "Individual", ylab= "Individual",main="Founder matches",
          labCol = geno_names, labRow = NA)
  dev.off()
}

