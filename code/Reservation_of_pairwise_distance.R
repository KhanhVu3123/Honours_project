source("FIt-SNE/fast_tsne.R",chdir = T)
# Load Fast-tSNE function

library(mvtnorm)
# to create multivariate normal distribution

library(Rtsne)
# to run normal tSNE.

sigma_6D <- rep(1,36)
i <- 1
while(i <= 36){
  sigma_6D[i] <- 2
  i <- i + 7
}
sigma_6D <- matrix(data = sigma_6D, nrow = 6, ncol = 6)

# Function to compute pairwise distance of the original data and dimensionality-reduced data and their correlation.
correlation_of_pairwise_distance <- function(original_df, tSNE_df){
  i <- 1
  
  # Calculate the pairwise distance for the original data.
  original_distance <- c()
  for(i in 1:(dim(original_df)[1] - 1)){
    num1 <- original_df[i,]
    for(j in i: dim(original_df)[1]){
      num2 <- original_df[j,]
      euc_dist <- sqrt(sum((num2 - num1)^2))
      original_distance <- c(original_distance, euc_dist)
    }
  }
  
  # Calculate the pairwise distance for the original data.
  tSNE_dist <- c()
  for(i in 1:(dim(tSNE_df)[1] - 1)){
    num1 <- tSNE_df[i,]
    for(j in i: dim(tSNE_df)[1]){
      num2 <- tSNE_df[j,]
      euc_dist <- sqrt(sum((num2 - num1)^2))
      tSNE_dist <- c(tSNE_dist, euc_dist)
    }
  }
  
  # Calculate Spearman and Pearson correlation coefficient.
  
  spearman_cor <- cor(original_distance, tSNE_dist, method = 'spearman')
  pearson_cor <- cor(original_distance, tSNE_dist, method = 'pearson')
  
  print(paste("Spearman correlation is", spearman_cor))
  print(paste("Pearson correlation is", pearson_cor))
}
# Initialize clusters in 6D to calculate the CPD.

dat6D_1 <- rmvnorm(n = 100, mean = rep(0,6), sigma = sigma_6D)
dat6D_2 <- rmvnorm(n = 100, mean = rep(10,6), sigma = sigma_6D)
dat6D_3 <- rmvnorm(n = 100, mean = rep(20,6), sigma = sigma_6D)

dat6D <- rbind(dat6D_1,dat6D_2,dat6D_3)

dat_tSNE <- Rtsne(dat6D, dims = 2)
correlation_of_pairwise_distance(dat6D, dat_tSNE$Y)

dat_fast_tSNE <- fftRtsne(dat6D, dims = 2)
correlation_of_pairwise_distance(dat6D, dat_fast_tSNE)

# Experiment was performed manunally for different number of clusters.

# Now plot them using ggplot.
library(ggplot2)

cluster_list <- c(rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4))
cor_method <- c(rep(c("tSNE_spearman","tSNE_pearson","fast_tSNE_spearman","fast_tSNE_pearson"),5))
cor_values <- c(0.897, 0.939, 0.861,0.913,0.764,0.797,0.851,0.936,0.516,0.579,0.926,0.981,0.551,0.580,0.731,0.778,0.560,0.570,0.693,0.749)

cor_df <- data.frame(Cluster = cluster_list, Value = cor_values, Method = cor_method)

pairwise_dist_plot <- ggplot(data = cor_df, aes(x = Cluster, y = Value,color = Method)) + geom_point() + ylim(0,1) + geom_line()

