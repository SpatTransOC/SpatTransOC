
library(tidyverse)
library(Seurat)
library(fgsea)
library(Scillus)

Idents(ov_combined) <- ov_combined$Group_Cluster

cluster_de <- list()

cluster_de[[1]] <- FindMarkers(ov_combined, ident.1 = 'PR_0', ident.2 = 'ER_0', 
                          logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)
cluster_de[[2]] <- FindMarkers(ov_combined, ident.1 = 'PR_1', ident.2 = 'ER_1', 
                               logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)
cluster_de[[3]] <- FindMarkers(ov_combined, ident.1 = 'PR_2', ident.2 = 'ER_2', 
                               logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)
cluster_de[[4]] <- FindMarkers(ov_combined, ident.1 = 'PR_3', ident.2 = 'ER_3', 
                               logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)
cluster_de[[5]] <- FindMarkers(ov_combined, ident.1 = 'PR_4', ident.2 = 'ER_4', 
                               logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)

Idents(ov_combined) <- ov_combined$Group

cluster_de[[6]] <- FindMarkers(ov_combined, ident.1 = 'PR', ident.2 = 'ER', 
                               logfc.threshold = -Inf, min.pct = 0, min.diff.pct = 0)


View(cluster_de[[1]])




de_to_stat(cluster_de[[1]], 'beta')

hist(de_to_stat(cluster_de[[2]], 'normal'))
hist(de_to_stat_bulk(res_bulk))

fgsea(pathways.hallmark, de_to_stat(cluster_de[[1]], 'beta'))
fgsea(pathways.hallmark, de_to_stat_bulk(res_bulk))

de_pathways <- list()
de_pathways[[1]] <- pathways.hallmark
de_pathways[[1]][['EMT-like']] <- kn_markers$`EMT-like`

pathways_tmp <- gmtPathways("~/Documents/projects/resources/c5.go.bp.v7.4.symbols.gmt")
de_pathways[[1]]$GOBP_CELL_POPULATION_PROLIFERATION <- pathways_tmp$GOBP_CELL_POPULATION_PROLIFERATION



cluster_de_gsea <- list()

for (i in 1:6) {
    cluster_de_gsea[[i]] <- fgsea(de_pathways[[1]], de_to_stat(cluster_de[[i]], mode = 'beta'))
}

cluster_de_gsea[[7]] <- fgsea(de_pathways[[1]], de_to_stat_bulk(res_bulk))

cluster_de_gsea <- do.call('rbind', cluster_de_gsea)

cluster_de_gsea$Comparison <- rep(c(paste0('PR_vs_ER_Cluster_', 
                                           c(as.character(c(0:4)), 'All')), 'PR_vs_ER_Bulk'), 
                                  each = length(de_pathways[[1]]))

# second set

pathway_names <- read_csv('./data/pathway_names.csv', col_names = F)[[1]]

pathways_tmp <- c(pathways_tmp,
                  gmtPathways("~/Documents/projects/resources/c2.all.v7.4.symbols.gmt"),
                  gmtPathways("~/Documents/projects/resources/c2.cp.reactome.v7.4.symbols.gmt"))

de_pathways[[2]] <- list()
for (i in pathway_names) {
    de_pathways[[2]][[i]] <- pathways_tmp[[i]]
}

cluster_de_gsea_2 <- list()

for (i in 1:6) {
    cluster_de_gsea_2[[i]] <- fgsea(de_pathways[[2]], de_to_stat(cluster_de[[i]], 'beta'))
}

cluster_de_gsea_2[[7]] <- fgsea(de_pathways[[2]], de_to_stat_bulk(res_bulk))

cluster_de_gsea_2 <- do.call('rbind', cluster_de_gsea_2)

cluster_de_gsea_2$Comparison <- rep(c(paste0('PR_vs_ER_Cluster_', 
                                           c(as.character(c(0:4)), 'All')), 'PR_vs_ER_Bulk'), 
                                  each = length(de_pathways[[2]]))



pdf("./figures/Cluster_group_GSEA_042121.pdf", width = 11, height = 11)

cluster_de_gsea %>%
    mutate(Significant = ifelse(padj < 0.1, 1, 0)) %>%
    ggplot(aes(x = pathway, y = Comparison)) +
    geom_point(aes(color = NES, size = -log10(padj))) +
    geom_point(aes(alpha = Significant), shape = 0, stroke = 1, size = 4, na.rm = T) +
    scale_alpha_continuous(range = c(0,1)) +
    scale_color_gradient2(low = '#3288bd', mid = 'white', high = '#d53e4f') +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(alpha = FALSE) +
    ggtitle('GSEA PR vs ER - 1')

cluster_de_gsea_2 %>%
    mutate(Significant = ifelse(padj < 0.1, 1, 0)) %>%
    ggplot(aes(x = pathway, y = Comparison)) +
    geom_point(aes(color = NES, size = -log10(padj))) +
    geom_point(aes(alpha = Significant), shape = 0, stroke = 1, size = 4, na.rm = T) +
    scale_alpha_continuous(range = c(0,1)) +
    scale_color_gradient2(low = '#3288bd', mid = 'white', high = '#d53e4f') +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(alpha = FALSE) +
    ggtitle('GSEA PR vs ER - 2')
dev.off()

hist(rbeta(1000,0.2,0.2)-0.5)

cluster_de_gsea %>%
    mutate(Significant = ifelse(padj < 0.1, 'Yes', ''))

# gsva

library(GSVA)

avg_exp_group <- AverageExpression(test, assays = 'SCT', group.by = 'group_cluster')[[1]]
clu_gsva <- gsva(avg_exp_group, pathways.hallmark)
cell_gsva <- gsva(as.matrix(exps), pathways.hallmark)

cell_gsva <- t(cell_gsva) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(group = test@meta.data$group,
           cluster = test@meta.data$seurat_clusters,
           group_cluster = test@meta.data$group_cluster)

cell_gsva <- cell_gsva %>%
    select(c(c(51:53), c(1:50)))

tmp <- cell_gsva %>%
    group_by(group_cluster) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)

as.matrix(tmp[2:51])[1:8,] - as.matrix(tmp[2:51])[10:17,]


# gene plots 052621

cancer_genes_1 <- read_csv("./data/cancer_census_list_1.csv", col_names = F)[[1]]

gene_de_df <- res_bulk[cancer_genes_1,] %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    mutate(Log2FC = log2FoldChange,
           PVal = pvalue,
           Comparison = 'PR_vs_ER_Bulk') %>%
    drop_na() %>%
    select(Gene, Log2FC, PVal, Comparison)


for (i in 0:4) {
    gene_de_df <- rbind(gene_de_df,
                        cluster_de[[i+1]][cancer_genes_1,] %>%
                            drop_na() %>%
                            rownames_to_column("Gene") %>%
                            mutate(Log2FC = avg_log2FC,
                                   PVal = p_val,
                                   Comparison = paste0('PR_vs_ER_Cluster_', i+1)) %>%
                            select(Gene, Log2FC, PVal, Comparison))
    
}

gene_de_df <- rbind(gene_de_df,
                    cluster_de[[6]][cancer_genes_1,] %>%
                        drop_na() %>%
                        rownames_to_column("Gene") %>%
                        mutate(Log2FC = avg_log2FC,
                               PVal = p_val,
                               Comparison = 'PR_vs_ER_Cluster_All') %>%
                        select(Gene, Log2FC, PVal, Comparison))

pdf("./figures/Gene_DE_052621.pdf", width = 8.5, height = 8.5)
gene_de_df %>%
    mutate(Significant = ifelse(PVal < 0.05, 1, 0)) %>%
    ggplot(aes(x = Comparison, y = Gene)) +
    geom_point(aes(size = -log10(PVal), color = ifelse(Log2FC > 0, 1, -1) * sqrt(abs(Log2FC)))) +
    geom_point(aes(alpha = Significant), shape = 0, stroke = 1, size = 4, na.rm = T) +
    scale_alpha_continuous(range = c(0,1)) +
    scale_color_gradient2(low = '#3288bd', mid = 'white', high = '#b2182b', midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(alpha = FALSE) +
    labs(color = 'Log2FC')
dev.off()


cancer_genes_2 <- read_csv("./data/cancer_census_list_2.csv", col_names = F)[[1]]


library(ComplexHeatmap)
library(circlize)

pdf("./figures/Gene_Cluster_Heatmap_052621.pdf", width = 8.5, height = 45)
ov_avg_exp_cluster[rownames(ov_avg_exp_cluster) %in% cancer_genes_2,] %>%
    Heatmap(col = colorRamp2(c(-1.5, 1, 3.5), c("#4575b4", "white", "#d73027")), name = "Expression")
dev.off()



