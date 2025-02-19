library(ggplot2)
source("FIt-SNE/fast_tsne.R",chdir = T)
library(umap)
library(cowplot)
library(patchwork)

head(fluorescent_downsampled)

# 1. Apply the standard arcsinh method with cofactor 150 to all the markers.
cofactor_set_1 <- c(5,150,250,500,1000)

cofactor <- cofactor_set_1[4]

fluorescent_downsampled_arcsinh <- asinh(fluorescent_downsampled/cofactor)
fluorescent_downsampled_umap <- umap(fluorescent_downsampled_arcsinh)

g1 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD3)) + geom_histogram()
g2 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD35)) + geom_histogram()
g3 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD4)) + geom_histogram()
g4 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD7)) + geom_histogram()
g5 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD8)) + geom_histogram()
g6 <- ggplot(data = fluorescent_downsampled_arcsinh, aes(CD27)) + geom_histogram()

plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3, ncol = 2)

arcsinh_150_grid <- list()
for(CDX in colnames(fluorescent_downsampled_arcsinh)){
  
  fluorescent_downsampled_umap_df <- data.frame(Dim1 = fluorescent_downsampled_umap$layout[,1],
                                                Dim2 = fluorescent_downsampled_umap$layout[,2],
                                                MarkerValue = fluorescent_downsampled_arcsinh[[CDX]])
  # We cannot make the CDX to be the column, instead put it to the title.
  g1 <- ggplot(fluorescent_downsampled_umap_df, aes(x = Dim1, y = Dim2, color = MarkerValue)) + geom_point() + scale_color_gradient(low = "blue", high = "red", name = CDX)
  arcsinh_150_grid[[CDX]] <- g1
}
arcsinh_150_plots <- wrap_plots(arcsinh_150_grid)
print(arcsinh_150_plots)

# This time we try with t-SNE.

fluorescent_downsampled_tsne <- fftRtsne(as.matrix(fluorescent_downsampled_arcsinh))
for(CDX in colnames(fluorescent_downsampled_arcsinh)){
  
  fluorescent_downsampled_tsne_df <- data.frame(Dim1 = fluorescent_downsampled_tsne[,1],
                                                Dim2 = fluorescent_downsampled_tsne[,2],
                                                MarkerValue = fluorescent_downsampled_arcsinh[[CDX]])
  # We cannot make the CDX to be the column, instead put it to the title.
  g1 <- ggplot(fluorescent_downsampled_tsne_df, aes(x = Dim1, y = Dim2, color = MarkerValue)) + geom_point() + scale_color_gradient(low = "blue", high = "red", name = CDX)
  arcsinh_150_grid[[CDX]] <- g1
}
arcsinh_150_plots <- wrap_plots(arcsinh_150_grid)
print(arcsinh_150_plots)

    
# CD35 and CD4 looks amazing!! Let's try other cofactor. We will try 5, 250, 500 and 1000 to see what happens.
# However, as I searched the marker up, I have multiple results for 1 marker.



