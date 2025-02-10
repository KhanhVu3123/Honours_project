library(cluster)
library(factoextra)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(ggpubr)
library(cowplot)

df <- read.csv("T_cell_fluorescent_size.csv")
fluorescent_signal <- df[,c(7,8,9,10,11,12)]
colnames(fluorescent_signal) <- c("CD3", "CD35", "CD4", "CD7", "CD8","CD27")

# Examine the histogram of all biomarkers in the original dataset
gg1 <- ggplot(data = fluorescent_signal, aes(x = CD3)) + geom_histogram(binwidth = 100) + xlim(c(-1000,1000))
gg2 <- ggplot(data = fluorescent_signal, aes(x = CD35)) + geom_histogram(binwidth = 100) + xlim(c(-100,10000))
gg3 <- ggplot(data = fluorescent_signal, aes(x = CD4)) + geom_histogram(binwidth = 100) + xlim(c(-100,20000))
gg4 <- ggplot(data = fluorescent_signal, aes(x = CD7)) + geom_histogram(binwidth = 100) + xlim(c(-100,20000))
gg5 <- ggplot(data = fluorescent_signal, aes(x = CD8)) + geom_histogram(binwidth = 100) + xlim(c(-500,5000))
gg6 <- ggplot(data = fluorescent_signal, aes(x = CD27)) + geom_histogram(binwidth = 100) + xlim(c(-300,1000))

plot_grid(gg1,gg2,gg3,gg4,gg5,gg6, labels = c("A","B","C","D","E","F"), nrow = 3, ncol = 2)

# Take 2000 random data points in the dataset and examine them. The aim is to see whether the downsampled-dataset still has the same properties as the original one.
random_index <- sample(1:(dim(fluorescent_signal)[1]), 2000)
fluorescent_signal_downsampled <- fluorescent_signal[random_index,]
fluorescent_signal_downsampled


for(colname in colnames(fluorescent_signal_downsampled)){
  
  g1 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD3)) + geom_histogram(binwidth = 100) + xlim(c(-1000,1000))
  g2 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD35)) + geom_histogram(binwidth = 100) + xlim(c(-100,10000))
  g3 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD4)) + geom_histogram(binwidth = 100) + xlim(c(-100,20000))
  g4 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD7)) + geom_histogram(binwidth = 100) + xlim(c(-100,20000))
  g5 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD8)) + geom_histogram(binwidth = 100) + xlim(c(-500,5000))
  g6 <- ggplot(data = fluorescent_signal_downsampled, aes(x = CD27)) + geom_histogram(binwidth = 100) + xlim(c(-300,1000))
  
  plot_grid(g1,g2,g3,g4,g5,g6, labels = c("A","B","C","D","E","F"), nrow = 3, ncol = 2)
  
}


