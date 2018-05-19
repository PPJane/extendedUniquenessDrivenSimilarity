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

# # 1. Goodall's similarity -----------------------------------------------
#### functions that return the Chi square transformed values
independent_nom_features_d <- function(x) {
  n <- length(x)
  freq <- sort(table(x))
  p2 <- freq*(freq-1)/(n*(n-1))
  
  cumfreq_freq <- c(0, cumsum(table(freq)))
  d <- cumsum(p2)
  for (i in 1:(length(cumfreq_freq)-1)) {
    d[(cumfreq_freq[i]+1):cumfreq_freq[i+1]] <- d[cumfreq_freq[i+1]]
  }
  
  equal_pairs <- names(d)
  names(d) <- paste(equal_pairs, equal_pairs)
  
  other_pairs <- combn(equal_pairs, 2)
  d_ext <- rep(1, ncol(other_pairs))
  names(d_ext) <- ifelse(other_pairs[1,] <= other_pairs[2,],
                         paste(other_pairs[1,], other_pairs[2,]),
                         paste(other_pairs[2,], other_pairs[1,]))
  
  missing_values <- paste(rep(0,length(equal_pairs)), equal_pairs)
  d_missing <- rep(1, length(missing_values))
  names(d_missing) <- missing_values
  
  d <- c(d, d_ext, d_missing)
  
  rep_counts <- table(d)
  n_unids <- length(rep_counts) # the number of unique d scores
  
  temp_d1 <- as.numeric(names(rep_counts)[-1]) # correspond to D_ij
  temp_d2 <- as.numeric(names(rep_counts)[-n_unids]) # correspond to (D_ij)'
  
  if (temp_d2[1]!=0){
    ki2 <- rep(c(2-2*log(temp_d2[1]), 2 - 2*(temp_d1*log(temp_d1) - temp_d2*log(temp_d2))/(temp_d1-temp_d2)), 
               rep_counts)
  }else{
    ki2 <- rep(c(Inf, c(2-2*log(temp_d1[2]), 2 - 2*(temp_d1[2:(n_unids-1)]*log(temp_d1[2:(n_unids-1)]) - temp_d2[2:(n_unids-1)]*log(temp_d2[2:(n_unids-1)]))/(temp_d1[2:(n_unids-1)]-temp_d2[2:(n_unids-1)]))), 
               rep_counts)
  }
  names(ki2) <- names(d)
  return(ki2)
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
  # a lot of tricks here
  # Firstly, for column diff_a_b, we compare the every pair of ith and (i+1)th value pairs
  # if they are equal, then it is possible that they are in the same group
  # We do further judgement by doing the same comparison for column uniqueness
  # the comparisons among the two columns will generate a m*2 logical matrix
  # By doing row sum, if the sum is 2, then the ith value pair and the (i+1)th value pair
  # should be in the same group, which means we need to use the last dist_a_b in this group
  # as the dist_a_b for other members
  # The part in the sort function gives an index vector where dist_a_b need to be modified
  # Say, if it is c(1, 2, 3, 5), then it means the 4th dist_a_b should be used for 
  # the 1st, 2nd, and 3rd member, the 6th dist_a_b should be used for the 5th member.
  # As a result, we can do a inverse loop and assigning values one by one
  # Note that this step can not be vectorized, because for the case of c(1,2,3), directly assigning
  # the 2nd dist_a_b to the 1st one is wrong, because the 2nd dist_a_b is also waiting to be 
  # corrected by itself. Here 2 should refer to 3 while 3 refer to 4. 
  # Thus we have to start from the end of the sequence and loop backwards
  rmat[,c(1,2)] <- c(univals[rmat[,1]], univals[rmat[,2]])
  # print(rmat)
  rownames(rmat) <- paste(rmat[,1], rmat[,2])
  
  ki2 <- -2*log(rmat[, 5])
  return(ki2) # return a named vector
}

dependent_nom_fts <- NULL
independent_nom_fts <- c('YelpTag1', 'YelpTag2', 'YelpTag3')
num_fts <- setdiff(colnames(mydata), c(dependent_nom_fts, independent_nom_fts))
t_c <- length(num_fts)
t_d <- length(independent_nom_fts) + length(dependent_nom_fts)
t_g <- 0

for (ft in colnames(mydata)){
  if (ft %in% num_fts){
    assign(paste0('ki2_', ft), num_features_d(mydata[[ft]]))
  } else {
    assign(paste0('ki2_', ft), independent_nom_features_d(mydata[[ft]]))
  }
}

for (ft in colnames(mydata)){
  assign(paste0('D_ki2_', ft), matrix(0, n, n))
}

for (ft in colnames(mydata)){
  temp <- get(paste0('D_ki2_', ft))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      ft_vals <- mydata[[ft]]
      val_a <- ft_vals[i]
      val_b <- ft_vals[j]
      val_pairs <- ifelse(val_a <= val_b,
                          paste(val_a, val_b),
                          paste(val_b, val_a))
      tmp <- get(paste0('ki2_', ft))
      temp[i,j] <- tmp[val_pairs]
    }
  }
  assign(paste0('D_ki2_', ft), temp)
}

Dlist_chi <- list(D_ki2_Amount, D_ki2_Rating, D_ki2_YelpTag1, D_ki2_YelpTag2, D_ki2_YelpTag3)