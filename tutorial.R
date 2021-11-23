
library(Seurat)
library(tidyverse)
library(patchwork)

e1 <- Load10X_Spatial("./data/II21476/outs/")
p1 <- Load10X_Spatial("./data/II21480/outs/")

e1$orig.ident <- "Extreme"
p1$orig.ident <- "Poor"

e1@images$slice1@key <- "E1_"
p1@images$slice1@key <- "P1_"

e1@images$E1 <- e1@images$slice1
e1@images$slice1 <- NULL

p1@images$P1 <- p1@images$slice1
p1@images$slice1 <- NULL

VlnPlot(e1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(e1, features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e1, features = "nCount_Spatial", pt.size.factor = 2.9) + theme(legend.position = "right")
SpatialFeaturePlot(e1, features = "nFeature_Spatial", pt.size.factor = 2.9) + theme(legend.position = "right")

VlnPlot(p1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(p1, features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p1, features = "nCount_Spatial", pt.size.factor = 2.8) + theme(legend.position = "right")
SpatialFeaturePlot(p1, features = "nFeature_Spatial", pt.size.factor = 2.8) + theme(legend.position = "right")

# subset

e1 <- subset(e1, slice1_imagecol <= 478)
p1 <- subset(p1, slice1_imagecol <= 415)


SpatialFeaturePlot(e1, features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right") 
SpatialFeaturePlot(e1, features = "nFeature_Spatial", pt.size.factor = 3) + theme(legend.position = "right")

SpatialFeaturePlot(p1, features = "nCount_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")
SpatialFeaturePlot(p1, features = "nFeature_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")

# Transform

e1 <- SCTransform(e1, assay = "Spatial", verbose = FALSE)
p1 <- SCTransform(p1, assay = "Spatial", verbose = FALSE)

# merge

ov <- merge(e1, p1)

DefaultAssay(ov) <- "SCT"
VariableFeatures(ov) <- c(VariableFeatures(e1), VariableFeatures(p1))
ov <- RunPCA(ov, verbose = FALSE)
ov <- FindNeighbors(ov, dims = 1:30)
ov <- FindClusters(ov, verbose = FALSE, resolution = 0.6)
ov <- RunUMAP(ov, dims = 1:30)


# Plot

wrap_plots(
    DimPlot(ov, reduction = "umap"),
    DimPlot(ov, reduction = "umap", group.by = 'orig.ident')
)

SpatialDimPlot(ov, label = TRUE, label.size = 3, pt.size.factor = 3)



# Markers


SpatialFeaturePlot(e1, 
                   features = head(SpatiallyVariableFeatures(e1, selection.method = "markvariogram"), 6), 
                   ncol = 3, alpha = c(0.1, 1), pt.size.factor = 3)

SpatialFeaturePlot(p1, 
                   str_subset(head(SpatiallyVariableFeatures(p1, selection.method = "markvariogram"), 13), "MT", negate = T), 
                   ncol = 3, alpha = c(0.1, 1), pt.size.factor = 3)



