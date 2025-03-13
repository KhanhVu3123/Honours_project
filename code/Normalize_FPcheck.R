# The aim in this script:

# 1. Checking whether normalizing the data will make any differerence to the clustering
# 2. Check for the false positive discovery of the FLowSOM tools.


library(ggplot2)
library(patchwork)
library(umap)
source("FIt-SNE/fast_tsne.R",chdir = T)
library(irlba)
library(cowplot)
library(patchwork)
library(uwot)
# For performing PCA. 

library(FlowSOM)
library(flowCore)
library(dplyr)

library(reshape2)
# To reshape data.

library(mvtnorm)
# To generate multivariate normal

library(mclust)

# Aim 1

head(fluorescent_signal_arcsinh)

# Try 2 different normalization technique, min-max and Z-score normalization

# Normalization 1: Min-max
# Function for normalizing a list
normalize <- function(df){
  return((df - min(df))/(max(df) - min(df)))
  
}

# Function for normalizing the whole df
normalize_df <- function(df){
  for(i in 1:ncol(df)){
    df[[i]] <- normalize(df[,i])
  }
  return(df)
}


test_df <- matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol =3)
test_df <- as.data.frame(test_df)

test_df_normalized <- normalize_df(test_df)
test_df_normalized

# Normalization 2: Z-score normalization
test_df_normalized_2 <- scale(test_df)
test_df_normalized_2

# Now perform on the real dataset

# Histogram plots of arcsinh-tranformed data

ggplot_list <- list()

for(CDx in colnames(fluorescent_signal_arcsinh)){
  g1 <- ggplot(fluorescent_signal_arcsinh, aes_string(x = CDx)) + geom_histogram()
  ggplot_list[[CDx]] <- g1
}
ggplot_wrap_plots <- wrap_plots(ggplot_list)
print(ggplot_wrap_plots)

## Min-max normalization 
fluorescent_signal_arcsinh_minmaxNorm <- normalize_df(fluorescent_signal_arcsinh)

# Visualize the output of min-max normalization
ggplot_list <- list()

for(CDx in colnames(fluorescent_signal_arcsinh_minmaxNorm)){
  g1 <- ggplot(fluorescent_signal_arcsinh_minmaxNorm, aes_string(x = CDx)) + geom_histogram()
  ggplot_list[[CDx]] <- g1
}
ggplot_wrap_plots <- wrap_plots(ggplot_list)
print(ggplot_wrap_plots)

## Z-score normalization
fluorescent_signal_arcsinh_Znorm <- scale(fluorescent_signal_arcsinh)

# Visualize the output of Z-score normalization
ggplot_list <- list()

for(CDx in colnames(fluorescent_signal_arcsinh_Znorm)){
  g1 <- ggplot(fluorescent_signal_arcsinh_Znorm, aes_string(x = CDx)) + geom_histogram()
  ggplot_list[[CDx]] <- g1
}
ggplot_wrap_plots <- wrap_plots(ggplot_list)
print(ggplot_wrap_plots)

# FlowSOM clustering testing

# Min-max testing
fSOM_minmax <- FlowSOM::ReadInput(as.matrix(fluorescent_signal_arcsinh_minmaxNorm), transform = F, scale = F)
fSOM_minmax <- FlowSOM::BuildSOM(fSOM_minmax, colsToUse = c(1,2,3,4,5,6))
fSOM_minmax <- FlowSOM::BuildMST(fSOM_minmax)

FlowSOM::PlotStars(fSOM_minmax)

# Z-score testing

fSOM_Znorm <- FlowSOM::ReadInput(as.matrix(fluorescent_signal_arcsinh_Znorm), transform = F, scale = F)
fSOM_Znorm <- FlowSOM::BuildSOM(fSOM_Znorm, colsToUse = c(1,2,3,4,5,6))
fSOM_Znorm <- FlowSOM::BuildMST(fSOM_Znorm)

FlowSOM::PlotStars(fSOM_Znorm)

# Min-max yields similar results, while in Z-score every cluster seems to express the same amount of every marker.
# Hypothesis: The negative value might cause the trouble, because the histogram of all markers are exactly the same !!


fluorescent_signal_arcsinh_Znorm_plus5 <- fluorescent_signal_arcsinh_Znorm + 5


fSOM_Znorm_plus5 <- FlowSOM::ReadInput(as.matrix(fluorescent_signal_arcsinh_Znorm_plus5), transform = F, scale = F)
fSOM_Znorm_plus5 <- FlowSOM::BuildSOM(fSOM_Znorm_plus5, colsToUse = c(1,2,3,4,5,6))
fSOM_Znorm_plus5 <- FlowSOM::BuildMST(fSOM_Znorm_plus5)

FlowSOM::PlotStars(fSOM_Znorm_plus5)

# Even when I add 5, making the majority of data positive, it still doesn't work.

# Aim 2: Test for false discovery of FlowSOM.

# Perform 2 tests: controlled synthetic data, and permutation testing

sigma_6D <- rep(1,36)
i <-1 
while(i <= 36){
  sigma_6D[i] <- 2
  i <- i + 7
}
sigma_6D <- matrix(sigma_6D, nrow = 6, ncol = 6)

# 2 clusters were twisted manually to investigate how well k-means does the clustering.
random_6Ddata <- rmvnorm(100000, mean = rep(5,6), sigma_6D)
random_6Ddata <- as.data.frame(random_6Ddata)


fSOM_random_6Ddata <- FlowSOM::ReadInput(as.matrix(random_6Ddata), transform = F, scale = F)
fSOM_random_6Ddata <- FlowSOM::BuildSOM(fSOM_random_6Ddata, colsToUse = c(1,2,3,4,5,6))
fSOM_random_6Ddata <- FlowSOM::BuildMST(fSOM_random_6Ddata)

FlowSOM::PlotStars(fSOM_random_6Ddata)

# kay, FlowSOM just spits out noise as expected.

# Now we shuffle things up.

shuffled_fluorescent_signal_arcsinh <- fluorescent_signal_arcsinh

for(i in 1:ncol(fluorescent_signal_arcsinh)){
  shuffled_fluorescent_signal_arcsinh[, i] <- sample(fluorescent_signal_arcsinh[,i])
}

# Visualize the shuffled data
ggplot_list <- list()

for(CDx in colnames(shuffled_fluorescent_signal_arcsinh)){
  g1 <- ggplot(shuffled_fluorescent_signal_arcsinh, aes_string(x = CDx)) + geom_histogram()
  ggplot_list[[CDx]] <- g1
}
ggplot_wrap_plots <- wrap_plots(ggplot_list)
print(ggplot_wrap_plots)

# The histogram of each marker still looks the same.

fSOM_shuffled_data <- FlowSOM::ReadInput(as.matrix(shuffled_fluorescent_signal_arcsinh), transform = F, scale = F)
fSOM_shuffled_data <- FlowSOM::BuildSOM(fSOM_shuffled_data, colsToUse = c(1,2,3,4,5,6))
fSOM_shuffled_data <- FlowSOM::BuildMST(fSOM_shuffled_data)

FlowSOM::PlotStars(fSOM_shuffled_data)

adjustedRandIndex(as.vector(fSOM_shuffled_data$map$mapping)[1:10000], as.vector(fSOM$map$mapping)[1:10000])
# Very close to 0, which is expected.

# However, the clusters in the shuffled dataset also makes perfect sense. If I didn't know that was shuffled data, I would still believe the results.



