source("FIt-SNE/fast_tsne.R",chdir = T)
library(umap)
library(mvtnorm)
library(hexbin)
library(RColorBrewer)
library(cluster)
library(scatterplot3d)
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


# Test with 1D

dat1D_1 <- rnorm(n = 1000, mean = 0 ,sd = 1)
dat1D_2 <- rnorm(n = 1000, mean = 5, sd = 2)

dat1D <- c(dat1D_1, dat1D_2)

hist(dat1D,breaks = 100, main = "")
plot(dat1D_df)


# Test with 2D
sigma2D_1 <- matrix(c(1,0.5,0.5,1), nrow = 2, ncol =2)
sigma2D_2 <- matrix(c(2,1,1,2), nrow = 2, ncol = 2)

dat2D_1 <- rmvnorm(n = 1000, mean = c(0,0), sigma = sigma2D_1)
dat2D_2 <- rmvnorm(n = 1000, mean = c(3,0), sigma = sigma2D_2)

dat2D <- rbind(dat2D_1, dat2D_2)
plot(dat2D, main = "Mu1 = (0,0), Mu2 = (3,0)")

kmeans_result_2D <- kmeans(dat2D, centers = 2)
gapstat <- clusGap(dat2D, FUN = kmeans, K.max = 10)
fviz_gap_stat(gapstat)
sum(kmeans_result_2D$cluster[1:100])
plot(dat2D, col = c(rep(1,100), rep(2,100)))
plot(kmeans_result_2D$cluster[1:100])

# Test with 3D
sigma3D_1 <- matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol =3)
sigma3D_2 <- matrix(c(3,0,0,0,3,0,0,0,3), nrow = 3, ncol = 3)


dat3D_1 <- rmvnorm(n = 1000, mean = c(0,0,0), sigma = sigma3D_1)
dat3D_2 <- rmvnorm(n = 1000, mean = c(3,0,0), sigma = sigma3D_2)

dat3D <- rbind(dat3D_1, dat3D_2)
scatterplot3d(dat3D)


dat3D_tSNE <- fftRtsne(dat3D, perplexity = 100)
plot(dat3D_tSNE)
kmeans_3D_ftSNE <- kmeans(dat3D_tSNE, centers = 2)
plot(dat3D_tSNE, col = kmeans_3D_ftSNE$cluster)
sum(kmeans_3D_ftSNE$cluster[seq(1,1000)])
# # Testing with 6D

sigma_6D <- rep(1,36)
i <-1 
while(i <= 36){
  sigma_6D[i] <- 2
  i <- i + 7
}
sigma_6D <- matrix(sigma_6D, nrow = 6, ncol = 6)

# 2 clusters were twisted manually to investigate how well k-means does the clustering.
dat6D_1 <- rmvnorm(100, mean = rep(5,6), sigma_6D)
dat6D_2 <- rmvnorm(100, mean = rep(9,6), sigma_6D)

dat6D <- rbind(dat6D_1, dat6D_2)
dat6D_std <- scale(dat6D)

dat6D_ftsne <- fftRtsne(dat6D)
plot(dat6D_ftsne)
# Clustering first 

# Using k-means
gapstat <- clusGap(dat6D_std, FUN = kmeans, K.max = 25)
print(fviz_gap_stat(gapstat))

dat6D_ftsne <- fftRtsne(dat6D, perplexity = 30)
print(plot(dat6D_ftsne, main = paste("Mu2 =",5)))


# Experiment with 3 clusters in 6D
acc_after <- c()
acc_before <- c()
calculate_mislabel <- function(clusters){
  most_appearances1 <- max(table(clusters[1001:2000]))
  most_appearances2 <- max(table(clusters[2001:3000]))
  most_appearances3 <- max(table(clusters[1:1000]))
  
  return(3000- most_appearances1 - most_appearances2 - most_appearances3)
}
for(run in 1:10){
  dat6D_1 <- rmvnorm(n = 1000, mean = rep(0 + run,6), sigma = sigma_6D)
  dat6D_2 <- rmvnorm(n = 1000, mean = rep(100 + run,6), sigma = sigma_6D)
  dat6D_3 <- rmvnorm(1000, mean = rep(103 + run,6), sigma = sigma_6D)
  #dat6D_4 <- rmvnorm(100, mean = rep(105,6), sigma = sigma_6D)
  
  
  dat6D <- rbind(dat6D_1, dat6D_2, dat6D_3)#, dat6D_4)
  # gapstat <- clusGap(scale(dat6D), FUN = kmeans, K.max = 10)
  # fviz_gap_stat(gapstat)
  kmeans_before_tSNE <- kmeans(dat6D, centers = 3)
  
  perplexity_set_1 <- c(1,10,30,50,60,70)
  par(mfrow = c(2,6))
  for (perplex in perplexity_set_1){
    dat6D_ftsne <- fftRtsne(dat6D)
    kmeans_after_tSNE <- kmeans(dat6D_ftsne, centers = 3)
    
  }
  
  acc1 <- calculate_mislabel(kmeans_after_tSNE$cluster)
  acc2 <- calculate_mislabel(kmeans_before_tSNE$cluster)
  
  acc_after <- c(acc_after, acc1)
  acc_before <- c(acc_before, acc2)
}



