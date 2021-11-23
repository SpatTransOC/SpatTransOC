
library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

for (i in 1:6) {
    e[[i]] <- SCTransform(e[[i]], assay = "Spatial", verbose = FALSE)
}

for (i in 1:5) {
    p[[i]] <- SCTransform(p[[i]], assay = "Spatial", verbose = FALSE)
}

for (i in 1:6) {
    e[[i]]$orig.ident <- paste0('E',i)
    e[[i]]@images$slice1@key <- paste0('E',i,'_')
}

for (i in 1:5) {
    p[[i]]$orig.ident <- paste0('P',i)
    p[[i]]@images$slice1@key <- paste0('P',i,'_')
}



e[[1]]@images$E1 <- e[[1]]@images$slice1
e[[2]]@images$E2 <- e[[2]]@images$slice1
e[[3]]@images$E3 <- e[[3]]@images$slice1
e[[4]]@images$E4 <- e[[4]]@images$slice1
e[[5]]@images$E5 <- e[[5]]@images$slice1
e[[6]]@images$E6 <- e[[6]]@images$slice1

e[[1]]@images$slice1 <- NULL
e[[2]]@images$slice1 <- NULL
e[[3]]@images$slice1 <- NULL
e[[4]]@images$slice1 <- NULL
e[[5]]@images$slice1 <- NULL
e[[6]]@images$slice1 <- NULL

p[[1]]@images$P1 <- p[[1]]@images$slice1
p[[2]]@images$P2 <- p[[2]]@images$slice1
p[[3]]@images$P3 <- p[[3]]@images$slice1
p[[4]]@images$P4 <- p[[4]]@images$slice1
p[[5]]@images$P5 <- p[[5]]@images$slice1


p[[1]]@images$slice1 <- NULL
p[[2]]@images$slice1 <- NULL
p[[3]]@images$slice1 <- NULL
p[[4]]@images$slice1 <- NULL
p[[5]]@images$slice1 <- NULL


# merge

test <- merge(e[[1]], e[[2]])

for (i in 3:6) {
    test <- merge(test, e[[i]])
}

for (i in 1:5) {
    test <- merge(test, p[[i]])
}

DefaultAssay(test) <- "SCT"
features <- unique(do.call("c", c(map(e, VariableFeatures), map(p, VariableFeatures))))


VariableFeatures(test) <- features
test <- RunPCA(test, verbose = FALSE)
test <- FindNeighbors(test, dims = 1:30)
test <- FindClusters(test, verbose = FALSE, resolution = 0.14)
test <- RunUMAP(test, dims = 1:30, min.dist = 0.05)

DimPlot(test, reduction = "umap", cols = c(brewer.pal(12, "Paired")))
DimPlot(test, reduction = "umap", group.by = 'orig.ident', cols = c(rep('blue',6), rep('red',5)))
DimPlot(test, reduction = "umap", group.by = 'orig.ident', cols = brewer.pal(12, "Paired"))

Idents(test) <- test@meta.data$seurat_clusters

pdf('./figures/Heatmap_041421.pdf', width = 11, height = 8.5)
DoHeatmap(test, do.call('c', map(.x = markers, .f = function(x) {x %>% head(10) %>% pull(Gene)})))
dev.off()

names(test@graphs)


pdf("./figures/Spatial_dim_plot_030721.pdf")
SpatialDimPlot(test, images = 'E1', cols = brewer.pal(12, 'Paired')[c(1:5,8)], pt.size.factor = 2)
SpatialDimPlot(test, images = 'E2', cols = brewer.pal(12, 'Paired')[1:8], pt.size.factor = 3)
SpatialDimPlot(test, images = 'E3', cols = brewer.pal(12, 'Paired')[c(1:4,6,8)], pt.size.factor = 4)
SpatialDimPlot(test, images = 'E4', cols = brewer.pal(12, 'Paired')[c(1:4,8)], pt.size.factor = 3.5)
SpatialDimPlot(test, images = 'E5', cols = brewer.pal(12, 'Paired')[c(1:4)], pt.size.factor = 2.5)
SpatialDimPlot(test, images = 'E6', cols = brewer.pal(12, 'Paired')[c(1:4,6,8)], pt.size.factor = 6)

SpatialDimPlot(test, images = 'P1', cols = brewer.pal(12, 'Paired')[c(1:4,7:9)], pt.size.factor = 3)
SpatialDimPlot(test, images = 'P2', cols = brewer.pal(12, 'Paired')[c(1,3,6,7)], pt.size.factor = 3)
SpatialDimPlot(test, images = 'P3', cols = brewer.pal(12, 'Paired')[c(1:4,6,8)], pt.size.factor = 3)
SpatialDimPlot(test, images = 'P4', cols = brewer.pal(12, 'Paired')[c(1:5)], pt.size.factor = 3)
SpatialDimPlot(test, images = 'P5', cols = brewer.pal(12, 'Paired')[c(1:4,6:8)], pt.size.factor = 3)
SpatialDimPlot(test, images = 'P6', cols = brewer.pal(12, 'Paired')[c(1:4,6,8)], pt.size.factor = 4)
dev.off()



