# install.packages('tidyr')
# install.packages('cocor')
# install.packages('pander')
library(dplyr)
library(tidyr)
library(cocor)
library(pander)

load('./data/cleanedSurveyData.rda')

likertAvg <- likertDF %>% group_by(SeedRestaurant, ComparisonRestaurant) %>% summarise(likertAvg = mean(LikertScore, na.rm = TRUE))
### some similarity measures below gives similarity scores above 1
### normalize all the similarity measures


distMeasurelist = list.files(path = './dissimilarityMatrices', pattern="*.csv")
seedrestaurantlist <- unique(likertAvg$SeedRestaurant)
seedComparison <- as.data.frame(matrix(0, nrow = 1+length(seedrestaurantlist), ncol = 2+length(distMeasurelist)*2))
colnames(seedComparison) <- c('Seed Restaurant', sapply(distMeasurelist, FUN = function(x){paste(x, c('estimate','p-value'), sep = '_')}), 'n')
seedComparison$`Seed Restaurant` <- c('All Seed Restaurants', seedrestaurantlist)

decimal_num <- 3
conf_level <- 0.95

k <- 1
for (i in distMeasurelist) {
  simScore <- algorithmDF %>% mutate(similarity = 1 - distance) %>% filter(distMeasure == i) %>%
    filter(distance != 0) %>% arrange(RestaurantName1, RestaurantName2)
    colnames(simScore) <- c('SeedRestaurant', 'ComparisonRestaurant', 'distMeasure', 'distance', 'similarity')
  comp <- merge(likertAvg, simScore, by = c('SeedRestaurant', 'ComparisonRestaurant'))
  
  tmp <- cor.test(comp$likertAvg, comp$similarity, alternative = 'two.sided', method = 'pearson', conf.level = conf_level)
  seedComparison[1, 2*k] <- round(tmp$estimate, digits = decimal_num)
  seedComparison[1, 2*k+1] <- round(tmp$p.value, digits = decimal_num)
  seedComparison[1, ncol(seedComparison)] <- tmp$parameter
  
  for (j in 1:length(seedrestaurantlist)) {
    tmp <- filter(comp, SeedRestaurant == seedrestaurantlist[j])
    tmp <- cor.test(tmp$likertAvg, tmp$similarity, alternative = 'two.sided', method = 'pearson', conf.level = conf_level)
    seedComparison[j+1, 2*k] <- round(tmp$estimate, digits = decimal_num)
    seedComparison[j+1, 2*k+1] <- round(tmp$p.value, digits = decimal_num)
    seedComparison[j+1, ncol(seedComparison)] <- tmp$parameter
  }
  k = k + 1
}



# visualize correlation significance --------------------------------------
for (i in grep(pattern = '_p-value', colnames(seedComparison))) {
  seedComparison[,i-1] <- paste(seedComparison[,i-1], add.significance.stars(seedComparison[,i]))
}
write.csv(seedComparison, file = 'seedComparison.csv')
