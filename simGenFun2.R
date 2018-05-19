# similarity generating functions
# 1. Goodall's similarity
# 2. Our proposed similarity
# 3. Our proposed distance similarity measure
# 4. Our proposed bag-of-tags similarity measure

library('dplyr')

# Read in data set --------------------------------------------------------
cat('\n*****Real Dataset*****\n')
mydata <- read.csv('./data/inputData.csv', stringsAsFactors = FALSE)
mydatalabel <- mydata$Restaurant
mydata <- mydata[,-1] # remove Restaurant Names
n <- nrow(mydata)
## numeric matrix feature Geographical Distance under walking mode
geodist <- read.csv("./data/address_walking.csv", stringsAsFactors = FALSE)
geodist <- geodist[,-1]
geodist[is.na(geodist)] <- 0
geodist <- geodist + t(geodist)

distance<-as.data.frame(table(unlist(geodist)))
distance[,3]<-distance$Freq/n/n
distance[,4] <- cumsum(distance[,3])
colnames(distance)=c("1Meter","2Freq","3Pr","4AccP")


# # 2. Our proposed similarity --------------------------------------------
dependent_nom_features_d <- function(x) {
  # x is going to be a matrix that cbind all the dependent features
  numfeatures <- ncol(x)
  n <- nrow(x)
  y <- as.vector(x)
  y <- y[-which(y=='0')]
  m <- length(y)
  freq <- sort(table(y))
  p2 <- freq*(freq-1)/(m*(m-1))
  
  cumfreq_freq <- c(0, cumsum(table(freq)))
  d <- cumsum(p2)
  for (i in 1:(length(cumfreq_freq)-1)) {
    d[(cumfreq_freq[i]+1):cumfreq_freq[i+1]] <- d[cumfreq_freq[i+1]]
  }
  
  equal_pairs <- names(d)
  names(d) <- paste(equal_pairs, equal_pairs)
  numuniqueval <- length(equal_pairs)
  
  instance_idx <- list()
  for (i in 1:numuniqueval){
    instance_idx[[i]] <- which(x==equal_pairs[i]) %% n
  }
  # find out for each unique value level, what instances include that value level
  
  co_occurrence <- NULL
  for (i in 1:(numuniqueval-1)){
    for (j in (i+1):numuniqueval){
      co_occurrence <- c(co_occurrence, length(intersect(instance_idx[[i]],instance_idx[[j]])))
    }
  }
  
  other_pairs <- combn(equal_pairs, 2)
  other_pairs <- other_pairs[,order(co_occurrence)]
  cumoccur_occur <- cumsum(table(co_occurrence))
  noncorrelated <- first(cumoccur_occur)
  other_pairs <- cbind(other_pairs[,-(1:noncorrelated)], 
                       other_pairs[,(1:noncorrelated)])
  cumoccur_occur <- c(cumoccur_occur[-1],cumoccur_occur[1])
  cumoccur_occur <- cumoccur_occur - noncorrelated
  cumoccur_occur[length(cumoccur_occur)] <- nth(cumoccur_occur, length(cumoccur_occur)-1) + noncorrelated
  cumoccur_occur <- c(0, cumoccur_occur)
  # reorder other_pairs based on co-occurrence, and move the noncorrelated pairs to the end of the matrix
  p2_ext <- 2*freq[other_pairs[1,]]*freq[other_pairs[2,]]/(m*(m-1))
  d_ext <- last(d) + cumsum(p2_ext)
  for (i in 1:(length(cumoccur_occur)-1)) {
    d_ext[(cumoccur_occur[i]+1):cumoccur_occur[i+1]] <- d_ext[cumoccur_occur[i+1]]
  }
  names(d_ext) <- ifelse(other_pairs[1,] <= other_pairs[2,],
                         paste(other_pairs[1,], other_pairs[2,]),
                         paste(other_pairs[2,], other_pairs[1,]))
  
  missing_values <- paste(rep(0,length(equal_pairs)+1), c('0',equal_pairs))
  d_missing <- rep(NA, length(missing_values))
  names(d_missing) <- missing_values
  
  d <- c(d, d_ext, d_missing)
  
  return(d)
}

num_features_d <- function(x) {
  freq <- table(x)
  univals <- as.numeric(dimnames(freq)[[1]]) 
  # unival stands for unique value, s for its plural form
  # note that table will return a frequency vector sorted by univals 
  # in the ascending order by default
  n_ob <- length(x)
  n <- length(univals)
  m <- n*(n+1)/2
  rmat <- matrix(nrow = m, ncol=5) # re-arranged matrix
  colnames(rmat) <- c('val_a', 'val_b', 'diff_a_b','uniqueness', 'dist_a_b')
  # the returned matrix
  # first assign 1:n
  rmat[1:n, c(1,2,4)] <- c(1:n, 1:n, freq)
  
  for (i in 1:(n-1)) {
    start_ind <- n*i - i*(i-1)/2 + 1 # the start index
    end_ind <- start_ind + n - i - 1 # the end index
    # print(start_ind, end_ind)
    pre_start_ind <- start_ind - (n-i+1) 
    pre_end_ind <- start_ind - 1
    rmat[start_ind:end_ind, c(1,2,4)] <- c(1:(n-i), 
                                           (i+1):n, 
                                           rmat[pre_start_ind:(pre_end_ind-1), 4] 
                                           + rmat[(i+1):n, 4])
  }
  rmat[,3] <- univals[rmat[,2]] - univals[rmat[,1]]
  rmat <- rmat[order(rmat[,3], rmat[,4]),]
  rmat[,5] <- cumsum(
    ifelse(
      rmat[,1] < rmat[,2], # val_a != val_b
      2*freq[rmat[,1]]*freq[rmat[,2]]/(n_ob*(n_ob-1)), # 2*f1*f2/(n*(n-1))
      freq[rmat[,1]]*(freq[rmat[,2]]-1)/(n_ob*(n_ob-1)) # f1*(f1-1)/(n*(n-1))
    )
  )
  # recall from nom_features_d, we need correct the results from naive cumsum
  to_corr_ind <- sort(
    (1:(m-1))[rowSums(rmat[1:(m-1),c(3,4)] == rmat[2:m, c(3,4)]) == 2], 
    decreasing = TRUE)
  for (ind in to_corr_ind) {
    rmat[ind, 5] <- rmat[ind + 1, 5]
  }
  rmat[,c(1,2)] <- c(univals[rmat[,1]], univals[rmat[,2]])
  rownames(rmat) <- paste(rmat[,1], rmat[,2])
  
  ki2 <- -2*log(rmat[, 5])
  return(ki2) # return a named vector
}


## numeric feature Geographical Distance by two addresses probability measure
dependent_nom_fts <- c('YelpTag1', 'YelpTag2', 'YelpTag3')
independent_nom_fts <- NULL
num_fts <- setdiff(colnames(mydata), c(dependent_nom_fts, independent_nom_fts))
t_c <- length(num_fts)
t_d <- length(independent_nom_fts) + 1 # length(dependent_nom_fts)
t_g <- 1

ki2_Amount <- num_features_d(mydata$Amount)
ki2_Rating <- num_features_d(mydata$Rating)
ki2_Tag <- dependent_nom_features_d(as.matrix(mydata[,dependent_nom_fts]))

D_ki2_Tag <- D_Tag <- D_ki2_Amount <- D_ki2_Rating <- D_ki2_Dist <- matrix(0, n, n)

for (i in 1:(n-1)){
  for (j in (i+1):n){
    D_ki2_Dist[i,j]=-2*log(distance[geodist[i,j]==distance$'1Meter',4])
  }
}

for (i in 1:(n-1)){
  for (j in (i+1):n){
    ft_vals <- mydata$Rating
    val_a <- ft_vals[i]
    val_b <- ft_vals[j]
    val_pairs <- ifelse(val_a <= val_b,
                        paste(val_a, val_b),
                        paste(val_b, val_a))
    D_ki2_Rating[i,j] <- ki2_Rating[val_pairs]
  }
}

for (i in 1:(n-1)){
  for (j in (i+1):n){
    val_a <- mydata$YelpTag1[i]
    val_b <- mydata[j,dependent_nom_fts]
    val_pairs <- ifelse(val_a <= val_b,
                        paste(val_a, val_b),
                        paste(val_b, val_a))
    tag1 <- ki2_Tag[val_pairs]
    
    val_a <- mydata$YelpTag2[i]
    val_b <- mydata[j,dependent_nom_fts]
    val_pairs <- ifelse(val_a <= val_b,
                        paste(val_a, val_b),
                        paste(val_b, val_a))
    tag2 <- ki2_Tag[val_pairs]
    
    val_a <- mydata$YelpTag3[i]
    val_b <- mydata[j,dependent_nom_fts]
    val_pairs <- ifelse(val_a <= val_b,
                        paste(val_a, val_b),
                        paste(val_b, val_a))
    tag3 <- ki2_Tag[val_pairs]
    
    tag_pair <- c(tag1, tag2, tag3)
    tag_pairs <- na.omit(tag_pair)
    D_Tag[i,j] <- weighted.mean(tag_pairs, w = exp(1-tag_pairs)/sum(exp(1-tag_pairs)))
  }
}

rep_counts <- table(D_Tag)
rep_counts <- as.numeric(names(rep_counts))
rep_counts <- c(0, rep_counts)

for (i in 1:(n-1)){
  for (j in (i+1):n){
    dp.idx <- which(abs(rep_counts - D_Tag[i,j]) < 0.0000001) - 1
    dp <- rep_counts[dp.idx]
    d <- D_Tag[i,j]
    if (d == 0){
      D_ki2_Tag[i,j] <- 2
    } else if (dp == 0){
      D_ki2_Tag[i,j] <- 2 - 2*log(d)
    } else {
      D_ki2_Tag[i,j] <- 2 - 2*((d*log(d)-dp*log(dp))/(d-dp))
    }
  }
}

Dlist_chi <- list(D_ki2_Tag, D_ki2_Amount, D_ki2_Rating, D_ki2_Dist)