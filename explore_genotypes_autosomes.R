#Load data
absdiff <- readRDS("Data/absdiff.rds")
attieDO <- readRDS("Data/base/attieDO_cross.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])

#generate column names
lcount <- 1
geno_names <- rep(NA,36)
letters <- c("A","B","C","D","E","F","G","H")
for(i in 1:8) {
  for(j in 1:8) {
    if(i <= j) {
      geno_names[lcount] <- paste0(letters[i],letters[j])
      lcount <- lcount + 1
    }
  }
}

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
        print(paste0(i,j))
        d[count] <- mean(gg == expected_geno(fgg[i,], fgg[j,]), na.rm=TRUE)
        count <- count + 1
      }
    }
  }
  
  d
}



##########################################################################################################
#Explore
library(viridis)

explore_markers <- c(6000,6403)
explore_ind <- seq(1,483)
explore_chr <- 3
round(matching_genos(ind=1,chr=explore_chr,markers=explore_markers),2)

for(i in explore_ind) {
  if(i==explore_ind[1]){
    matching <- round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2)
  }  
  else {
    matching <- rbind(matching, round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2))
  }
}

colnames(matching) <- geno_names
matching

explore_ind_names <- ind_names[explore_ind]
heatmap(matching,Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,1),
        xlab = "Individual", ylab= "Individual",main="Founder matches",
        labCol = geno_names)






explore_markers <- c(200,1200)
explore_ind <- seq(1,483)
explore_chr <- 16
round(matching_genos(ind=explore_ind,chr=explore_chr,markers=explore_markers),2)

for(i in explore_ind) {
  if(i==explore_ind[1]){
    matching <- round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2)
  }  
  else {
    matching <- rbind(matching, round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2))
  }
}

colnames(matching) <- geno_names
matching

explore_ind_names <- ind_names[explore_ind]
heatmap(matching,Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,1),
        xlab = "Individual", ylab= "Individual",main="Founder matches",
        labCol = geno_names)






#all markers
explore_ind <- seq(1,483)
explore_chr <- 3
explore_markers <- c(1,dim(absdiff[[explore_chr]])[2])
round(matching_genos(ind=explore_ind,chr=explore_chr,markers=explore_markers),2)

for(i in explore_ind) {
  if(i==explore_ind[1]){
    matching <- round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2)
  }  
  else {
    matching <- rbind(matching, round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2))
  }
}

colnames(matching) <- geno_names
matching

explore_ind_names <- ind_names[explore_ind]
heatmap(matching,Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,1),
        xlab = "Individual", ylab= "Individual",main="Founder matches",
        labCol = geno_names)








explore_markers <- c(200,800) #compare to 200-800
explore_ind <- seq(1,242)
explore_chr <- 20

for(i in explore_ind) {
  if(i==explore_ind[1]){
    matching <- round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2)
  }  
  else {
    matching <- rbind(matching, round(matching_genos(ind=i,chr=explore_chr,markers=explore_markers),2))
  }
}

colnames(matching) <- geno_names
matching

explore_ind_names <- ind_names[explore_ind]
heatmap(matching,Colv=NA, Rowv=NA, col = viridis(256), zlim=c(0,1),
        xlab = "Individual", ylab= "Individual",main="Founder matches",
        labCol = geno_names)


