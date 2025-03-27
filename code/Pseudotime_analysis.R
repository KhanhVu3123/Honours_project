
fluorescent_signal_arcsinh_matrix <- as.matrix(fluorescent_signal_arcsinh)
cell_ind <- seq(1,372385,1)
cell_ind <- paste0("c",cell_ind)

rownames(fluorescent_signal_arcsinh_matrix) <- cell_ind

fluorescent_signal_arcsinh_matrix_transposed <- t(fluorescent_signal_arcsinh_matrix)

sce_myDat <- SingleCellExperiment(assays = List(counts = fluorescent_signal_arcsinh_matrix_transposed))
reducedDims(sce_myDat) <- SimpleList(tSNE = fluorescent_signal_tSNE, UMAP = fluorescent_signal_umap)

clusterLabel_matrix <- as.matrix(fluorescent_signal_arcsinh_with_labels$Cluster)
rownames(clusterLabel_matrix) <- cell_ind

colData(sce_myDat)$SOM <- clusterLabel_matrix

sce_myDat_results <- slingshot(sce_myDat, clusterLabels = 'SOM', reducedDim = 'UMAP')

sce_myDat_results_tSNE <- slingshot(sce_myDat, clusterLabels = 'SOM', reducedDim = 'tSNE')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_myDat_results$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_myDat_results)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_myDat_results), lwd=2, col = 'black')#,type = 'lineages', col='black')


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_myDat_results_tSNE$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_myDat_results)$tSNE, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_myDat_results), lwd=2, col = 'black')#,type = 'lineages', col='black')






labels <- fluorescent_signal_with_labels_RNAEx
colors <- brewer.pal(9, 'Set1')[labels]

plot(reducedDims(sce_myDat_results)$UMAP, 
     col = colors[match(sce_myDat_results$SOM, labels)], 
     pch = 16, asp = 1)

# Add legend
legend("topright",                # or "bottomright", etc.
       legend = fluorescent_signal_with_labels_RNAEx$CellType,          # text labels
       col = brewer.pal(9, 'Set1')[labels],  # matching colors
       pch = 16,
       title = "SOM Cluster")

lines(SlingshotDataSet(sce_myDat_results), lwd=2, col = 'black')

lineages <- getLineages(data = dimred, clusterLabels = clustering)

