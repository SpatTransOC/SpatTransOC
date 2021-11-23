
library(GSVA)
library(ComplexHeatmap)
library(circlize)

ov_combined$Group_Cluster <- paste0(ov_combined$Group, "_", ov_combined$seurat_clusters)

table(ov_combined$Group_Cluster)

Idents(ov_combined)

div_pal <- colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027"))

ov_avg_exp_cluster <- AverageExpression(ov_combined, assays = 'integrated')[[1]]

Idents(ov_combined) <- ov_combined$Group_Cluster

ov_avg_exp_group_cluster <- AverageExpression(ov_combined, assays = 'integrated')[[1]]

# for LR analysis
ov_avg_exp_group_cluster <- AverageExpression(ov_combined, assays = 'SCT')[[1]]


pdf('./figures/Cluster_GSVA_042121.pdf')
Heatmap(mtx_percentile(gsva(ov_avg_exp_cluster, kn_markers)), name = "Score", col = div_pal)

Heatmap(mtx_percentile(gsva(ov_avg_exp_group_cluster[,str_detect(colnames(ov_avg_exp_group_cluster), 'E')], 
                            kn_markers)),
        col = div_pal,
        name = 'Score',
        cluster_columns = F, column_title = 'ER') +
Heatmap(mtx_percentile(gsva(ov_avg_exp_group_cluster[,str_detect(colnames(ov_avg_exp_group_cluster), 'P')], 
                                kn_markers)),
        col = div_pal,
        name = 'Score',
        cluster_columns = F, column_title = 'PR')
dev.off()


# barplot 052621

library(GSVA)

sample(names(ov_combined$orig.ident), size = 1000)

gsva_list <- list()

Idents(ov_combined) <- ov_combined$Group_Cluster

for (i in 1:50) {
        set.seed(i)
        cells <- sample(names(ov_combined$orig.ident), size = 1000)
        tmp <- subset(ov_combined, cells = cells)
        gsva_list[[i]] <- gsva(AverageExpression(tmp, assays = 'integrated', 
                                                 slot = 'scale.data')[[1]],
                                        kn_markers)
}

rm(cells, tmp)

gsva_list[[1]]
gsva_list[[2]]


gsva_reshape <- function(mtx) {
        
        mtx <- mtx %>% 
                as.data.frame() %>%
                rownames_to_column('Pathway') %>%
                pivot_longer(!Pathway, names_to = 'Group_Cluster', values_to = 'Score')
        
        return(mtx)
}

gsva_df <- purrr::map(gsva_list, gsva_reshape)

gsva_df <- do.call('rbind', gsva_df)

gsva_df <- gsva_df %>%
        mutate(Group = str_extract(Group_Cluster, '.*R'),
               Cluster = str_extract(Group_Cluster, '[0-9]'))

library(ggpubr)

pdf('./figures/Cluster_GSVA_052621.pdf', width = 44, height = 33)
gsva_df %>%
        filter(Cluster %in% c(0:5,8)) %>%
        arrange(Cluster) %>%
        ggboxplot(fill = 'Group', y = 'Score', x = 'Cluster', palette = 'Set2') +
        stat_compare_means(aes(group = Group), label = "p.format") +
        facet_wrap(~Pathway, scales = 'free') 
dev.off()


               