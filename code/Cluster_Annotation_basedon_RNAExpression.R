library(tidyverse)

# Aim: Annotate the clusters based on RNA expression from ImmGen microarray 

# 1. Annotate clusters by calculating the minimum Euclidean distance from the clusters to the "ground truth"


RNA_expression <- read.csv("RNA_expression_abTcells.csv")

head(RNA_expression)

par(mfrow = c(3,2))

for(col_ind in 2:ncol(RNA_expression)){
  hist(RNA_expression[[col_ind]], main = colnames(RNA_expression)[col_ind])
}

for(col_ind in 2:ncol(RNA_expression)){
  RNA_expression[[col_ind]] <- asinh(RNA_expression[[col_ind]]/150)
}

RNA_expression <- RNA_expression[1:49,]
min_dist_vector <- c()

CellType_label <- c()
for(i in 1:nrow(cluster_expression_df)){
  expression_vector <- as.numeric(cluster_expression_df[i, 2:6])
  distance_vector <- c()
  for(j in 1:nrow(RNA_expression)){
    cell_type_expression_vector <- as.numeric(RNA_expression[j, 2:6])
    euc_dist <- sqrt(sum((expression_vector - cell_type_expression_vector)^2))
    distance_vector <- c(distance_vector, euc_dist)
  }
  min_ind <- which(distance_vector == min(distance_vector))
  min_dist_vector <- c(min_dist_vector, min(distance_vector))
  CellType_label <- c(CellType_label, RNA_expression$CellType[min_ind])
}

cluster_expression_df$CellType <- CellType_label

# 2. UMAP/tSNE plot of the clusters

annotation <- cluster_expression_df[,c(1,8)]

annotation$Labels <- as.numeric(annotation$Labels)
annotation <- annotation %>% arrange(annotation$Labels)
annotation$CellType[11] <- "Dead"
annotation$CellType[29] <- "Dead"

fluorescent_signal_with_labels_RNAEx <- fluorescent_signal_arcsinh_with_labels[,1:7]
colnames(fluorescent_signal_with_labels_RNAEx)[colnames(fluorescent_signal_with_labels_RNAEx) == "Labels"] <- "Cluster"

fluorescent_signal_with_labels_RNAEx$Cluster <- as.numeric(fluorescent_signal_with_labels_RNAEx$Cluster)
fluorescent_signal_with_labels_RNAEx <- fluorescent_signal_with_labels_RNAEx %>% left_join(annotation, by = c("Cluster" = "Labels"))

fluorescent_signal_tSNE_to_plot$Labels <- fluorescent_signal_with_labels_RNAEx$CellType[rand_ind]
colnames(fluorescent_signal_tSNE_to_plot) <- c("tSNE_X", "tSNE_Y", "Labels")
ggplot(fluorescent_signal_tSNE_to_plot, aes(x = tSNE_X, y = tSNE_Y, color = as.factor(Labels))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                                         "#FFFF33", "#A65628", "#00CED1", "#999999", "#66C2A5",
                                         "#FC8D62", "#8DA0CB", "#E7298A", "#A6D854")) +
                                           theme_minimal()

# Explanations of cell type labels:
# 1. NKT.4-Lv: Natural killer cells that are CD4- from the liver.
# 2. NKT.4-Sp: Natural killer cells that are CD4- from the spleen
# 3. NKT.44-NK1.1-Th: Thymus pre-NKT CD44- and NK1.1-
# 4. T.4+8int.Th: DP that prepares to be CD8 cells.
# 5. T.4FP3+25+.Sp: Tregs (regulatory T cells) that are CD25+ from the spleen
# 6. T.4SP24int.Th: Semi-mature CD4+ T cells
# 7. T.DP69+.Th: DP, early positive selection.
# 





