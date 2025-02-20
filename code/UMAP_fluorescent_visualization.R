library(ggplot2)
source("FIt-SNE/fast_tsne.R",chdir = T)
library(umap)
library(cowplot)
library(patchwork)

head(fluorescent_downsampled)
# Examine the random downsampled data again (this has 2000 data points out of 372000)

binwidth <- 100

g1 <- ggplot(data = fluorescent_downsampled, aes(CD3)) + geom_histogram(binwidth = 10) + xlim(range(-100,2000))
g2 <- ggplot(data = fluorescent_downsampled, aes(CD35)) + geom_histogram(binwidth = 20) + xlim(range(-100,5000))
g3 <- ggplot(data = fluorescent_downsampled, aes(CD4)) + geom_histogram(binwidth = binwidth) + xlim(range(-100,15000))
g4 <- ggplot(data = fluorescent_downsampled, aes(CD7)) + geom_histogram(binwidth = binwidth)+ xlim(range(-500,10000))
g5 <- ggplot(data = fluorescent_downsampled, aes(CD8)) + geom_histogram(binwidth = binwidth) + xlim(range(-100,2000))
g6 <- ggplot(data = fluorescent_downsampled, aes(CD27)) + geom_histogram(binwidth = 10) + xlim(range(-100,300))

plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3, ncol = 2)




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

# 2. Apply the linlog transforamation method 

linlog_theta_set <- c(1,5,10,20)

linlog_theta <- linlog_theta_set[4]
fluorescent_downsampled_linlog <- linlog_for_df(fluorescent_downsampled, linlog_theta)

# Visualize the transformed data
binwidth <- 0.5
g1 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD3)) + geom_histogram(binwidth = binwidth) + xlim(range(0,12))
g2 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD35)) + geom_histogram(binwidth = binwidth)+ xlim(range(-100,15))
g3 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD4)) + geom_histogram(binwidth = binwidth)+ xlim(range(0,10))
g4 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD7)) + geom_histogram(binwidth = binwidth)+ xlim(range(0,10))
g5 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD8)) + geom_histogram(binwidth = binwidth)+ xlim(range(0,10))
g6 <- ggplot(data = fluorescent_downsampled_linlog, aes(CD27)) + geom_histogram(binwidth = binwidth)+ xlim(range(-50,10))

plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3, ncol = 2)
warnings()
# Perform UMAP, tSNE and visualize 
fluorescent_downsampled_umap <- umap(fluorescent_downsampled_linlog)
fluorescent_downsampled_tsne <- fftRtsne(as.matrix(fluorescent_downsampled_linlog))

fluorescent_dimension_reduction <- fluorescent_downsampled_tsne
linlog_1_grid <- list()
for(CDX in colnames(fluorescent_downsampled_linlog)){
  if(length(class(fluorescent_dimension_reduction)) == 1){
    fluorescent_downsampled_umap_df <- data.frame(Dim1 = fluorescent_dimension_reduction$layout[,1],
                                                  Dim2 = fluorescent_dimension_reduction$layout[,2],
                                                  MarkerValue = fluorescent_downsampled[[CDX]])
  }
  else{
    fluorescent_downsampled_umap_df <- data.frame(Dim1 = fluorescent_dimension_reduction[,1],
                                                  Dim2 = fluorescent_dimension_reduction[,2],
                                                  MarkerValue = fluorescent_downsampled[[CDX]])
    
  }
  # We cannot make the CDX to be the column, instead put it to the title.
  g1 <- ggplot(fluorescent_downsampled_umap_df, aes(x = Dim1, y = Dim2, color = MarkerValue)) + geom_point() + scale_color_gradient(low = "blue", high = "red", name = CDX)
  linlog_1_grid[[CDX]] <- g1
}
linlog_1_grid <- wrap_plots(linlog_1_grid)
print(linlog_1_grid)

# 1 question arises: Can we use one parameter for this column and another parameter for different column !?

# 3.Apply the Box-Cox transformation

boxcox_theta_set <- c(-5,-3,-1,1,3)

boxcox_theta <- boxcox_theta_set[3]
fluorescent_downsampled_boxcox <- boxcox_transformation_for_df(fluorescent_downsampled, boxcox_theta)

g1 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD3)) + geom_histogram(binwidth = binwidth)# + xlim(range(0,12))
g2 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD35)) + geom_histogram(binwidth = binwidth)# + xlim(range(-100,15))
g3 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD4)) + geom_histogram(binwidth = binwidth)# + xlim(range(0,10))
g4 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD7)) + geom_histogram(binwidth = binwidth)# + xlim(range(0,10))
g5 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD8)) + geom_histogram(binwidth = binwidth)# + xlim(range(0,10))
g6 <- ggplot(data = fluorescent_downsampled_boxcox, aes(CD27)) + geom_histogram(binwidth = binwidth)# + xlim(range(-50,10))

plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3, ncol = 2)
warnings()

# Perform UMAP, tSNE and visualize 
fluorescent_downsampled_umap <- umap(fluorescent_downsampled_boxcox)
fluorescent_downsampled_tsne <- fftRtsne(as.matrix(fluorescent_downsampled_boxcox))

fluorescent_dimension_reduction <- fluorescent_downsampled_tsne
boxcox_1_grid <- list()
for(CDX in colnames(fluorescent_downsampled_boxcox)){
  if(length(class(fluorescent_dimension_reduction)) == 1){
    fluorescent_downsampled_umap_df <- data.frame(Dim1 = fluorescent_dimension_reduction$layout[,1],
                                                  Dim2 = fluorescent_dimension_reduction$layout[,2],
                                                  MarkerValue = fluorescent_downsampled[[CDX]])
  }
  else{
    fluorescent_downsampled_umap_df <- data.frame(Dim1 = fluorescent_dimension_reduction[,1],
                                                  Dim2 = fluorescent_dimension_reduction[,2],
                                                  MarkerValue = fluorescent_downsampled[[CDX]])
    
  }
  # We cannot make the CDX to be the column, instead put it to the title.
  g1 <- ggplot(fluorescent_downsampled_umap_df, aes(x = Dim1, y = Dim2, color = MarkerValue)) + geom_point() + scale_color_gradient(low = "blue", high = "red", name = CDX)
  boxcox_1_grid[[CDX]] <- g1
}
boxcox_1_grid <- wrap_plots(boxcox_1_grid)
print(boxcox_1_grid)




