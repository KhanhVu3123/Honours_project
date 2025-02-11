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
  sigma_6D[i] <- 2
  i <- i + 7
}
sigma_6D <- matrix(sigma_6D, nrow = 6, ncol = 6)

# 2 clusters were twisted manually to investigate how well k-means does the clustering.
dat6D_1 <- rmvnorm(100, mean = rep(5,6), sigma_6D)
dat6D_2 <- rmvnorm(100, mean = rep(7,6), sigma_6D)

dat6D <- rbind(dat6D_1, dat6D_2)
dat6D_std <- scale(dat6D)

# Clustering first 

# Using k-means
gapstat <- clusGap(dat6D_std, FUN = kmeans, K.max = 25)
print(fviz_gap_stat(gapstat))

dat6D_ftsne <- fftRtsne(dat6D, perplexity = 30)
print(plot(dat6D_ftsne, main = paste("Mu2 =",5)))



# Experiment with 3 clusters in 6D

dat6D_1 <- rmvnorm(n = 100, mean = rep(0,6), sigma = sigma_6D)
dat6D_2 <- rmvnorm(n = 100, mean = rep(100,6), sigma = sigma_6D)
dat6D_3 <- rmvnorm(100, mean = rep(103,6), sigma = sigma_6D)
dat6D_4 <- rmvnorm(100, mean = rep(105,6), sigma = sigma_6D)


dat6D <- rbind(dat6D_1, dat6D_2, dat6D_3, dat6D_4)
gapstat <- clusGap(scale(dat6D), FUN = kmeans, K.max = 10)
fviz_gap_stat(gapstat)
kmeans_before_tSNE <- kmeans(dat6D, centers = 4)

perplexity_set_1 <- c(1,10,30,50,60,70)
par(mfrow = c(2,6))
for (perplex in perplexity_set_1){
  dat6D_ftsne <- fftRtsne(dat6D)
  kmeans_after_tSNE <- kmeans(dat6D_ftsne, centers = 4)
  plot(dat6D_ftsne, main = paste0("Perplexity = ",perplex,", before tSNE"), col = kmeans_before_tSNE$cluster)
  plot(dat6D_ftsne, main = paste0("Perplexity = ",perplex, ", after tSNE"), col = kmeans_after_tSNE$cluster)
  
}


