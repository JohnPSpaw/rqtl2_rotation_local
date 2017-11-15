jump_intervals <- function(chr,ind,threshold) {
  #matrix for chromosome
  matrix <- absdiff[[chr]];
  
  if(exists("endpoints")) {remove(endpoints)}
  #probs are abs diff at each marker ---->  set all probs below threshold to 0
  probs <- matrix[ind,]; probs[probs < threshold] <- 0; #print(probs)
  markers <- seq(1,nmar) 
  
  #marker numbers at which jumps occur
  jump_markers <- c(0,which(abs(diff(probs))>=threshold),length(probs)) ; #print(jump_markers)
  
  #lengths of intervals (alternating between below and above threshold)
  lengths <- diff(jump_markers); #print(lengths)
  
  #Plot thresholded probs with red lines at jump points
  plot(probs, main=paste0("Ind ",ind)) 
  abline(v=jump_markers, col="red")
  
  #create list of endpoints
  if(probs[1]==0){   #first prob below threshold
    indices <- seq(2,length(jump_markers),2)
    for(i in indices) {
      left <- jump_markers[i];  right <- jump_markers[i+1]
      if(i==2) {endpoints <- c(left,right)} #condition for first element to create dataframe
      else {endpoints <- rbind(endpoints,c(left,right))}
    }
  }
  else {
    indices <- seq(1,length(jump_markers),2)
    for(i in indices) {
      left <- jump_markers[i];  right <- jump_markers[i+1]
      if(i==1) {endpoints <- c(left,right)}
      else {endpoints <- rbind(endpoints,c(left,right))}
    }
  }
  
  print(endpoints)
  #Only continue with individuals with meaningful jump points
  condition <- !((length(endpoints)==2)&is.na(endpoints[2])) 
  if(condition) {
    #Change column names in array
    #dimnames(endpoints)[[2]] <- c("Left", "Right")
    return(endpoints)    
  }
}

for(k in 60:41) {
  writeLines(paste0("\n","Ind ",k))
  jump_intervals(chr=10,ind=k,threshold=1.75)
}













