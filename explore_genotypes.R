library(viridis)
#Load data
absdiff <- readRDS("Data/absdiff.rds")
attieDO <- readRDS("Data/base/attieDO_cross.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])



explore_markers <- c(6000,6403)
explore_ind <- seq(1,483)
explore_chr <- 3
explore_genos(explore_chr,explore_ind,explore_markers)

explore_ind <- seq(1,483)
explore_chr <- 3
explore_markers <- c(1,dim(absdiff[[explore_chr]])[2])
explore_genos(explore_chr,explore_ind,explore_markers)





explore_markers <- c(4100,4600)
explore_ind <- seq(1,483)
explore_chr <- 8
explore_genos(explore_chr,explore_ind,explore_markers)

explore_ind <- seq(1,483)
explore_chr <- 8
explore_markers <- c(1,dim(absdiff[[explore_chr]])[2])
explore_genos(explore_chr,explore_ind,explore_markers)





explore_markers <- c(200,1200)
explore_ind <- seq(1,483)
explore_chr <- 16
explore_genos(explore_chr,explore_ind,explore_markers)

explore_ind <- seq(1,483)
explore_chr <- 16
explore_markers <- c(1,dim(absdiff[[explore_chr]])[2])
explore_genos(explore_chr,explore_ind,explore_markers)







explore_markers <- c(200,1000) #compare to 200-800
explore_ind <- seq(1,242)
explore_chr <- 20
explore_genos(explore_chr,explore_ind,explore_markers)

explore_ind <- seq(1,242)
explore_chr <- 20
explore_markers <- c(1,dim(absdiff[[explore_chr]])[2])
explore_genos(explore_chr,explore_ind,explore_markers)



