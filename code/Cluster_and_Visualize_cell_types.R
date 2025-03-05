# Plan for this script:
# 1. Logicle transformation
# 2. Arcsinh transformation with different cofactor
# 3. Most importantly: FlowSOM clustering


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

head(fluorescent_signal)
colnames(fluorescent_signal) <- c("CD44", "CD8", "CD25", "CD28", "CD4", "live_dead_viability")

# Do tranformation with cofactor 150.
set.seed(30102003)

fluorescent_signal_arcsinh <- asinh(fluorescent_signal/150)
fluorescent_signal_arcsinh_5 <- asinh(fluorescent_signal/5)
fluorescent_signal_logicle <- flowCore::biexponentialTransform(as.matrix(fluorescent_signal))

# Do flowSOM clustering

fSOM <- FlowSOM::ReadInput(as.matrix(fluorescent_signal_arcsinh), transform = F, scale = F)
fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = c(1,2,3,4,5,6))
fSOM <- FlowSOM::BuildMST(fSOM)

FlowSOM::PlotStars(fSOM)

str(fSOM$map)
head(fSOM$map$mapping)
dim(fSOM$map$mapping)


# Perform metaclustering, by hand, we know that we have 10 different cell types here (DN1, DN3a, DN3b, DN3c, DN4a, DN4b, DP, CD4 and CD8)
metaClustering <- as.character(metaClustering_consensus(fSOM$map$codes[c(-1,-91,-96),], k =10))
metaClustering <- append(metaClustering, "Unk", after = 0)
metaClustering <- append(metaClustering, "Unk", after = 90)
metaClustering <- append(metaClustering, "Unk", after = 95)

PlotLabels(fSOM, labels = metaClustering)

# Label the data from clustering to metaclusters(cell types)
labels <- fSOM$map$mapping[,1]
clusters <- c()

for(label in labels){
  clusters <- c(clusters, metaClustering[label])
}

count <- 0

cell_type_labels <- c()
for(cluster in clusters){
  if(count %% 1000 == 0){
    print(count)
  }
  
  if(cluster == 1){
    cell_type_labels <- c(cell_type_labels,"DN3")
  }
  if(cluster %in% c(2,3,4)){
    cell_type_labels <- c(cell_type_labels, "CD8")
  }
  if(cluster == 5){
    cell_type_labels <- c(cell_type_labels, "DN1/2")
  }
  if(cluster == 6){
    cell_type_labels <- c(cell_type_labels, "DeadCells")
  }
  if(cluster == 7){
    cell_type_labels <- c(cell_type_labels, "CD4")
  }
  if(cluster == 8){
    cell_type_labels <- c(cell_type_labels, "DN4")
  }
  if(cluster %in% c(9,10)){
    cell_type_labels <- c(cell_type_labels, "DP")
  }
  if(cluster == "Unk"){
    cell_type_labels <- c(cell_type_labels, "Unk")
  }
  count <- count + 1
  
}

# Visualize the clusters in 2D using UMAP and tSNE 

# tSNE 
fluorescent_signal_tSNE <- fftRtsne(as.matrix(fluorescent_signal_arcsinh), max_iter = 3000, learning_rate = as.integer(300000/12))


fluorescent_signal_tSNE_to_plot <- fluorescent_signal_tSNE[rand_ind,]
fluorescent_signal_tSNE_to_plot <- as.data.frame(fluorescent_signal_tSNE_to_plot)

colnames(fluorescent_signal_tSNE_to_plot) <- c("tSNE_X", "tSNE_Y")
fluorescent_signal_tSNE_to_plot$labels <- as.vector(cell_type_labels[rand_ind])
colnames(fluorescent_signal_tSNE_to_plot) <- c("tSNE_X", "tSNE_Y", "Labels")
ggplot(fluorescent_signal_tSNE_to_plot, aes(x = tSNE_X, y = tSNE_Y, color = Labels)) + geom_point()


# UMAP
fluorescent_signal_umap <- uwot::umap(fluorescent_signal_arcsinh)
fluorescent_signal_umap_to_plot <- fluorescent_signal_umap[rand_ind,]
fluorescent_signal_umap_to_plot <- as.data.frame(fluorescent_signal_umap_to_plot)
colnames(fluorescent_signal_umap_to_plot) <- c("UMAP_X", "UMAP_Y")
fluorescent_signal_umap_to_plot$Labels <- as.vector(cell_type_labels[rand_ind])
ggplot(fluorescent_signal_umap_to_plot, aes(x = UMAP_X, y = UMAP_Y, color = Labels)) + geom_point()

# Plotting heatmap of expression for each cluster.

fluorescent_signal_arcsinh_with_labels <- fluorescent_signal_arcsinh
fluorescent_signal_arcsinh_with_labels$Labels <- cell_type_labels

cell_type_expression_df <- data.frame()

for(cell_type in unique(cell_type_labels)){
  unique_cell_type <- fluorescent_signal_arcsinh_with_labels %>% filter(Labels == cell_type)
  
  temp_df <- data.frame(Labels = cell_type)
  
  for(colname in colnames(unique_cell_type)[1:6]){
    temp_df[[colname]] <- mean(unique_cell_type[[colname]])
  }
  cell_type_expression_df <- bind_rows(cell_type_expression_df, temp_df)
  print(cell_type_expression_df)
}

heatmap_data <- melt(cell_type_expression_df, id.vars = "Labels")
ggplot(data= heatmap_data, aes(x = variable, y = Labels, fill = value)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) 

# Visualize the clusters without the unknown and dead cells
fluorescent_signal_arcsinh_without_unk_deadcells <- fluorescent_signal_arcsinh_with_labels %>% filter(Labels != "Unk" & Labels != "DeadCells")
fluorescent_signal_umap <- uwot::umap(fluorescent_signal_arcsinh_without_unk_deadcells)
rand_ind <- sample(1:nrow(fluorescent_signal_umap), 20000)
fluorescent_signal_umap_to_plot <- fluorescent_signal_umap[rand_ind,]
fluorescent_signal_umap_to_plot <- as.data.frame(fluorescent_signal_umap_to_plot)
colnames(fluorescent_signal_umap_to_plot) <- c("UMAP_X", "UMAP_Y")

cell_type_labels <- as.vector(cell_type_labels)
cell_type_labels_without_unk_deadcells <- cell_type_labels

for (type_to_remove in c("DeadCells","Unk")){
  cell_type_labels_without_unk_deadcells <- cell_type_labels_without_unk_deadcells[cell_type_labels_without_unk_deadcells != type_to_remove]
  
}
fluorescent_signal_umap_to_plot$Labels <- as.vector(cell_type_labels_without_unk_deadcells[rand_ind])
ggplot(fluorescent_signal_umap_to_plot, aes(x = UMAP_X, y = UMAP_Y, color = Labels)) + geom_point()




