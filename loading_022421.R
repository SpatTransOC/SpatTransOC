
library(Seurat)
library(tidyverse)
library(patchwork)

e <- list()
p <- list()

e[[1]] <- Load10X_Spatial("./data/II21472/outs/")
e[[2]] <- Load10X_Spatial("./data/II21473/outs/")
e[[3]] <- Load10X_Spatial("./data/II21474/outs/")
e[[4]] <- Load10X_Spatial("./data/II21475/outs/")
e[[5]] <- Load10X_Spatial("./data/II21476/outs/")
e[[6]] <- Load10X_Spatial("./data/II21477/outs/")

p[[1]] <- Load10X_Spatial("./data/II21478/outs/")
p[[2]] <- Load10X_Spatial("./data/II21480/outs/")
p[[3]] <- Load10X_Spatial("./data/II21481/outs/")
p[[4]] <- Load10X_Spatial("./data/II21482/outs/")
p[[5]] <- Load10X_Spatial("./data/II21483/outs/")

# remove outliers

e[[1]] <- subset(e[[1]], cells = setdiff(rownames(e[[1]]@meta.data), outliers[[1]]))
e[[2]] <- subset(e[[2]], cells = setdiff(rownames(e[[2]]@meta.data), outliers[[2]]))
e[[3]] <- subset(e[[3]], cells = setdiff(rownames(e[[3]]@meta.data), outliers[[3]]))
e[[4]] <- subset(e[[4]], cells = setdiff(rownames(e[[4]]@meta.data), outliers[[4]]))
e[[5]] <- subset(e[[5]], cells = setdiff(rownames(e[[5]]@meta.data), outliers[[5]]))
e[[6]] <- subset(e[[6]], cells = setdiff(rownames(e[[6]]@meta.data), outliers[[6]]))

p[[1]] <- subset(p[[1]], cells = setdiff(rownames(p[[1]]@meta.data), outliers[[7]]))
p[[2]] <- subset(p[[2]], cells = setdiff(rownames(p[[2]]@meta.data), outliers[[8]]))
p[[3]] <- subset(p[[3]], cells = setdiff(rownames(p[[3]]@meta.data), outliers[[9]]))
p[[4]] <- subset(p[[4]], cells = setdiff(rownames(p[[4]]@meta.data), outliers[[10]]))
p[[5]] <- subset(p[[5]], cells = setdiff(rownames(p[[5]]@meta.data), outliers[[11]]))






SpatialFeaturePlot(e[[1]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[1]], features = "nCount_Spatial", pt.size.factor = 2.7) + theme(legend.position = "right")
SpatialFeaturePlot(e[[1]], features = "nFeature_Spatial", pt.size.factor = 2.7) + theme(legend.position = "right")

SpatialFeaturePlot(e[[2]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[2]], features = "nCount_Spatial", pt.size.factor = 3.7) + theme(legend.position = "right")
SpatialFeaturePlot(e[[2]], features = "nFeature_Spatial", pt.size.factor = 3.7) + theme(legend.position = "right")

SpatialFeaturePlot(e[[3]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[3]], features = "nCount_Spatial", pt.size.factor = 5) + theme(legend.position = "right")
SpatialFeaturePlot(e[[3]], features = "nFeature_Spatial", pt.size.factor = 5) + theme(legend.position = "right")

SpatialFeaturePlot(e[[4]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[4]], features = "nCount_Spatial", pt.size.factor = 4) + theme(legend.position = "right")
SpatialFeaturePlot(e[[4]], features = "nFeature_Spatial", pt.size.factor = 4) + theme(legend.position = "right")

SpatialFeaturePlot(e[[5]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[5]], features = "nCount_Spatial", pt.size.factor = 2.9) + theme(legend.position = "right")
SpatialFeaturePlot(e[[5]], features = "nFeature_Spatial", pt.size.factor = 2.9) + theme(legend.position = "right")

SpatialFeaturePlot(e[[6]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(e[[6]], features = "nCount_Spatial", pt.size.factor = 8) + theme(legend.position = "right")
SpatialFeaturePlot(e[[6]], features = "nFeature_Spatial", pt.size.factor = 8) + theme(legend.position = "right")

SpatialFeaturePlot(p[[1]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[1]], features = "nCount_Spatial", pt.size.factor = 3.6) + theme(legend.position = "right")
SpatialFeaturePlot(p[[1]], features = "nFeature_Spatial", pt.size.factor = 3.6) + theme(legend.position = "right")

SpatialFeaturePlot(p[[2]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[2]], features = "nCount_Spatial", pt.size.factor = 4) + theme(legend.position = "right")
SpatialFeaturePlot(p[[2]], features = "nFeature_Spatial", pt.size.factor = 4) + theme(legend.position = "right")

SpatialFeaturePlot(p[[3]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[3]], features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right")
SpatialFeaturePlot(p[[3]], features = "nFeature_Spatial", pt.size.factor = 3) + theme(legend.position = "right")

SpatialFeaturePlot(p[[4]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[4]], features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right")
SpatialFeaturePlot(p[[4]], features = "nFeature_Spatial", pt.size.factor = 3) + theme(legend.position = "right")

SpatialFeaturePlot(p[[5]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[5]], features = "nCount_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")
SpatialFeaturePlot(p[[5]], features = "nFeature_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")

SpatialFeaturePlot(p[[6]], features = "nCount_Spatial", pt.size.factor = 0) + theme(legend.position = "right")
SpatialFeaturePlot(p[[6]], features = "nCount_Spatial", pt.size.factor = 5) + theme(legend.position = "right")
SpatialFeaturePlot(p[[6]], features = "nFeature_Spatial", pt.size.factor = 5) + theme(legend.position = "right")


# process








