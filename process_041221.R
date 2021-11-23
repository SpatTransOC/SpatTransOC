
library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(Scillus)

ov <- list()

ov[[1]] <- Load10X_Spatial("./data/II21472/outs/", slice = 'ER_1')
ov[[2]] <- Load10X_Spatial("./data/II21473/outs/", slice = 'ER_2')
ov[[3]] <- Load10X_Spatial("./data/II21474/outs/", slice = 'ER_3')
ov[[4]] <- Load10X_Spatial("./data/II21475/outs/", slice = 'ER_4')
ov[[5]] <- Load10X_Spatial("./data/II21476/outs/", slice = 'ER_5')
ov[[6]] <- Load10X_Spatial("./data/II21477/outs/", slice = 'ER_6')

ov[[7]] <- Load10X_Spatial("./data/II21478/outs/", slice = 'PR_1')
ov[[8]] <- Load10X_Spatial("./data/II21480/outs/", slice = 'PR_2')
ov[[9]] <- Load10X_Spatial("./data/II21481/outs/", slice = 'PR_3')
ov[[10]] <- Load10X_Spatial("./data/II21482/outs/", slice = 'PR_4')
ov[[11]] <- Load10X_Spatial("./data/II21483/outs/", slice = 'PR_5')

for (i in 1:11) {
    ov[[i]]$Sample <- c(paste0('ER_',1:6),paste0('PR_',1:5))[i]
    ov[[i]]$Group <- c(rep('ER',6),rep('PR',5))[i]
}


# remove outliers

for (i in 1:11) {
    ov[[i]] <- subset(ov[[i]], cells = setdiff(rownames(ov[[i]]@meta.data), outliers[[i]]))
}

# CCA

ov <- lapply(X = ov, FUN = SCTransform, assay = 'Spatial')

ov_features <- SelectIntegrationFeatures(ov, nfeatures = 15000)

ov <- PrepSCTIntegration(object.list = ov, anchor.features = ov_features)

saveRDS(ov, "./rdata/ov.rds")
saveRDS(ov_features, "./rdata/ov_features.rds")

ov_combined <- IntegrateData(FindIntegrationAnchors(ov, 
                                                    normalization.method = "SCT", 
                                                    anchor.features = ov_features,
                                                    assay = rep('SCT',11)),
                             normalization.method = "SCT")

saveRDS(ov_combined, "./rdata/ov_combined_041221.rds")

ov_combined <- readRDS("./data/ov_combined.rds")

ov_combined <- RunPCA(ov_combined)
ElbowPlot(ov_combined, ndims = 50)
ov_combined <- RunUMAP(ov_combined, 
                       reduction = "pca", 
                       dims = 1:50,
                       n.neighbors = 70)

ov_combined <- FindNeighbors(ov_combined, reduction = "pca", dims = 1:30)
ov_combined <- FindClusters(ov_combined, resolution = 0.4)

ov_combined@meta.data

plot_scdata(ov_combined, pal_setup = c("Paired"))
plot_scdata(ov_combined, pal_setup = "Set2", color_by = "Group")
plot_scdata(ov_combined, pal_setup = c("Set2","Set1"), color_by = "Sample")

pdf("./figures/Cluster_Stat_042121.pdf")
plot_stat(ov_combined, 'group_count', group_by = 'seurat_clusters', pal_setup = "Paired")
plot_stat(ov_combined, 'prop_fill', group_by = 'Group', pal_setup = 'Paired')
plot_stat(ov_combined, 'prop_fill', group_by = 'Sample', pal_setup = 'Paired')
dev.off()

ov_markers <- FindAllMarkers(ov_combined, logfc.threshold = 0.5, min.pct = 0.1, only.pos = T)

pdf("./figures/Heatmap_052621.pdf", width = 13, height = 8.5)
plot_heatmap(ov_combined, 
             get_top_genes(ov_combined,ov_markers,8)[1:48],
             n = 8, 
             anno_var = 'seurat_clusters', 
             anno_colors = list('Paired'),
             hm_limit = c(-1,0,1))
dev.off()



# RPCA

ov <- lapply(X = ov, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

ov_rpca_features <- SelectIntegrationFeatures(object.list = ov)

ov <- lapply(X = ov, FUN = function(x) {
    x <- ScaleData(x, features = ov_rpca_features, verbose = FALSE)
    x <- RunPCA(x, features = ov_rpca_features, verbose = FALSE)
})

ov_rpca <- IntegrateData(FindIntegrationAnchors(object.list = ov, 
                                                anchor.features = ov_rpca_features,
                                                reduction = "rpca"))

ov_rpca <- ScaleData(ov_rpca, verbose = FALSE)
ov_rpca <- RunPCA(ov_rpca, npcs = 30, verbose = FALSE)
ov_rpca <- RunUMAP(ov_rpca, reduction = "pca", dims = 1:30)
ov_rpca <- FindNeighbors(ov_rpca, reduction = "umap", dims = 1:2)
ov_rpca <- FindClusters(ov_rpca, resolution = 0.2)

ov_rpca@meta.data

DimPlot(ov_rpca, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(ov_rpca, reduction = "umap", group.by = 'Sample', label = TRUE)
DimPlot(ov_rpca, reduction = "umap", group.by = 'Group', label = TRUE)
DimPlot(ov_rpca, reduction = "umap", label = TRUE, split.by = 'Group', ncol = 2)

plot_stat(ov_rpca, 'prop_fill', group_by = 'Group')

ov_rpca_markers <- FindAllMarkers(ov_rpca, logfc.threshold = 0.2, min.pct = 0, only.pos = T)

DoHeatmap(ov_rpca, ov_rpca_markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC) %>% pull(gene))






