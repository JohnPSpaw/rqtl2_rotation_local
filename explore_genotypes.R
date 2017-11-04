################################################################
#LOAD DATA
#Cross object
attieDO <- readRDS("Data/base/attieDO_cross.rds")
# raw genotypes
g <-attieDO$geno[["X"]] # matrix individuals x markers
# founder genotypes
fg <-attieDO$founder_geno[["X"]] # matrix founders x markers
# first 200 markers
markers <- colnames(g)[1:1000]
# sexes
is_female <- attieDO$is_female
females <- which(is_female)
males <- which(!is_female)

################################################################
matching_genos <- function(ind, chr, markers) {
  g <-attieDO$geno[[chr]] 
  fg <-attieDO$founder_geno[[chr]]

  #markers
  left_endpoint <- markers[1]; right_endpoint <- markers[2];
  markers <- colnames(g)[left_endpoint:right_endpoint]
  #print(length(markers))
  #genotypes
  gg <- g[ind, markers]
  fgg <- fg[,markers]
  # problem: missing values are coded as 0
  gg[gg==0] <- NA
  fgg[fgg==0] <- NA
  
  #print(apply(fgg, 1, function(a) mean(a == gg, na.rm=TRUE)))
}


#Manually examine genotype proportions at specified markers for specified individuals
inds <- c(58,75)
markers <- c(1500,2200)
for(k in inds[1]:inds[2]){ 
  #print(paste0("Ind ",k))
  if(k==inds[1]){genos <- matching_genos(ind=k, chr="X", markers)}
  else {genos <- rbind(genos, matching_genos(ind=k, chr="X", markers))}
}
rownames(genos) <- paste0(seq(inds[1],inds[2]))
print(genos)

#Manually examine genotype proportions at specified markers for specified individuals
inds <- c(30,60)
markers <- c(250,700)
for(k in inds[1]:inds[2]){ 
  #print(paste0("Ind ",k))
  if(k==inds[1]){genos <- matching_genos(ind=k, chr="X", markers)}
  else {genos <- rbind(genos, matching_genos(ind=k, chr="X", markers))}
}
rownames(genos) <- paste0(seq(inds[1],inds[2]))
print(genos)