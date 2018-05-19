# chi square transformation -----------------------------------------------
chisqAgg <- function(Dlist){
  for (i in Dlist){
    if (isSymmetric(i)){
      i[lower.tri(i)] <- 0
    }
  }
  D_ki2 <- Reduce("+", Dlist) / 2
  D <- matrix(0,n,n)
  k_seq <- 0:(t_d + t_c + t_g - 1) # used for \sum_{k=0}^{t_c+t_d-1}
  fac_seq <- gamma(k_seq+1) # factorial(x) (x! for non-negative integer x)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      D[i,j] <- exp(-D_ki2[i,j])*sum(D_ki2[i,j]**k_seq/fac_seq)
    }
  }
  D <- D + t(D)
  return(D)
}

D <- chisqAgg(Dlist_chi)

# save the aggregated dissimilarity matrx ---------------------------------
# run the Goodall's similarity and our proposed similarity measure separately;
# then save the dissimilarity matrices separately by uncomment the corresponding line and run it.
# # 1. Goodall's similarity measure
# write.csv(D, file = './dissimilarityMatrices/simFun1_chisqAgg.csv', row.names = mydatalabel)
# # 2. Our proposed similarity measure
# write.csv(D, file = './dissimilarityMatrices/simFun2_chisqAgg.csv', row.names = mydatalabel)
# # 3. Our proposed similarity measure - Bag-of-tags
# write.csv(D, file = './dissimilarityMatrices/simFun3_chisqAgg.csv', row.names = mydatalabel)
# # 4. Our proposed similarity measure - Distance Matrix
# write.csv(D, file = './dissimilarityMatrices/simFun4_chisqAgg.csv', row.names = mydatalabel)