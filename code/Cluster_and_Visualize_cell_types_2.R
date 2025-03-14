# In the first version of this script, I performed meta-clustering with 8 cell types in total, not accounting for the intermediate cell types (e.g. ISP CD8+)
# Here, I will increase the number of clusters in the clusterin part and not allowing any more unknown in my dataset.

# Aim: Annotate the clusters with 30 meta-cluster.


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

PlotStars(fSOM)

metaClustering <- as.character(metaClustering_consensus(fSOM$map$codes, k =30))
save(fSOM, file = "initial_clustering_results.RData")

PlotLabels(fSOM, labels = metaClustering)

labels <- fSOM$map$mapping[,1]
clusters <- c()

for(label in labels){
  clusters <- c(clusters, metaClustering[label])
}

fluorescent_signal_arcsinh_with_labels <- fluorescent_signal_arcsinh
fluorescent_signal_arcsinh_with_labels$Labels <- clusters


cluster_expression_df <- data.frame()

for(cluster in unique(clusters)){
  cluster_population <- fluorescent_signal_arcsinh_with_labels %>% filter(Labels == cluster)
  
  temp_df <- data.frame(Labels = cluster)
  
  for(colname in colnames(cluster_population)[1:6]){
    temp_df[[colname]] <- mean(cluster_population[[colname]])
  }
  cluster_expression_df <- bind_rows(cluster_expression_df, temp_df)
  print(cluster_expression_df)
}

cluster_expression_df$Labels <- factor(cluster_expression_df$Labels, levels = dplyr::mixedsort(unique(cluster_expression_df$Labels))) # Order numerically


heatmap_data <- melt(cluster_expression_df, id.vars = "Labels")
heatmap_data$Labels <- factor(heatmap_data$Labels, levels = seq(1,30))
ggplot(data= heatmap_data, aes(x = variable, y = Labels, fill = value)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) 



