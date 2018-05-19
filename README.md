# extendedUniquenessDrivenSimilarity

code and data release for paper: A Uniqueness-Driven Similarity Measure for Automated Competitor Identification
scripts wrote in R
Data is released under data subfolder in csv and RData files.

simGenFun1-4.R each corresponding to one similarity calculation:
  1. Goodall's similarity
  2. Our proposed similarity measure
  3. Our proposed distance similarity measure
  4. Our proposed bag-of-tags similarity measure

Instruction to replicate: 
  1. Set your working directory to where you saved the R scripts and data
  2. Source the four files separately
  3. For each source, run simAggregation.R (chisqAgg function to aggregate the dissimilarity calculated for each feature), and then save the similarity matrix by uncomment and run the corresponding write.csv command.
  4. Once you have generated the dissimilarity matrix file, source likertComparison.R file to check the results, a seedComparison.csv will be generated in your working directory.
