source("FIt-SNE/fast_tsne.R",chdir = T)
library(umap)
library(mvtnorm)

library(cluster)
# to perform k-means clustering

library(factoextra)
# to perform visualization of k-means results.

library(fpc)
# Experiment design:
# 1. How many dimensions? A: 6
# 2. How many clusters? A: We test it with 2 to 6 like the last time
# 3. How far should each cluster be? Also test it with different distances.
# 4. What algorithms should I use? Should be unsupervised clusering e.g. Louvain, clustering, k-means, KNN??
# 5. How many data points do we need? Total of 400,000.
# 6. What accuracy metrics to be use here? A: CPD, KNN!?

# Testing k-means

sigma_6D <- rep(2,36)
i <-1 
while(i <= 36){
  sigma_6D[i] <- 4
  i <- i + 7
}
sigma_6D <- matrix(sigma_6D, nrow = 6, ncol = 6)

# 2 clusters were twisted manually to investigate how well k-means does the clustering.
dat6D_1 <- rmvnorm(100, mean = rep(100,6), sigma_6D)
dat6D_2 <- rmvnorm(100, mean = rep(104,6), sigma_6D)

dat6D <- rbind(dat6D_1, dat6D_2)
dat6D_std <- scale(dat6D)

# Clustering first 

# Using k-means
gapstat <- clusGap(dat6D_std, FUN = kmeans, K.max = 25)
fviz_gap_stat(gapstat)



