

# Markers

ov_markers <- list()

for (i in 0:5) {
    ov_markers[[i+1]] <- FindMarkers(ov_combined, ident.1 = i, logfc.threshold = 0.5)
}

View(markers[[1]])
View(markers[[2]])
View(markers[[3]])
View(markers[[4]])
View(markers[[5]])
View(markers[[6]])
View(markers[[7]])
View(markers[[8]])
View(markers[[9]])

for (i in 1:6) {
    ov_markers[[i]] <- ov_markers[[i]] %>%
        rownames_to_column('Gene') %>%
        arrange(-avg_log2FC)
}

do.call('c', map(.x = ov_markers, .f = function(x) {x %>% head(10) %>% pull(Gene)}))

names(markers) <- paste0("Cluster_",0:8)

library(openxlsx)

write.xlsx(markers, file = "./tables/markers_030721.xlsx")

SpatialFeaturePlot(e1, 
                   features = rownames(head(markers[[1]] %>% filter(avg_log2FC > 0), 6)), 
                   ncol = 3, alpha = 1, pt.size.factor = 3)

SpatialFeaturePlot(e1, 
                   features = rownames(head(markers[[2]] %>% filter(avg_log2FC > 0), 6)), 
                   ncol = 3, alpha = 1, pt.size.factor = 3)

SpatialFeaturePlot(e1, 
                   features = rownames(head(markers[[3]] %>% filter(avg_log2FC > 0), 6)), 
                   ncol = 3, alpha = 1, pt.size.factor = 3)


# all markers

avg_exp <- AverageExpression(test, assays = 'SCT')[[1]]
colnames(avg_exp) <- paste0('Cluster_', colnames(avg_exp))

avg_exp %>%
    as.data.frame() %>%
    rownames_to_column('Gene') %>%
    as_tibble() %>%
    write_tsv(file = './tables/all_markers_030921.tsv')

avg_exp %>%
    as.data.frame() %>%
    rownames_to_column('Gene') %>%
    as_tibble() %>%
    pivot_longer(cols = !Gene, 
                 names_to = 'Cluster', 
                 values_to = 'Expression') %>%
    arrange(desc(Expression)) %>%
    group_by(Cluster) %>%
    slice(1:5000) %>%
    write_tsv(file = './tables/top25pct_markers_030921.tsv')


# known marker scoring

mtx_percentile <- function(mtx) {
    
    mtx2 <- mtx
    
    for (i in 1:ncol(mtx)) {
        mtx2[,i] <- (rank(mtx[,i])/(length(mtx[,i])+1) - 0.5) * 2
    }
    return(mtx2)
}

kn_markers <- read_csv('./data/cell_type_makers_040721.csv')
kn_markers <- kn_markers %>%
    split(.,.['Cell type'])
kn_markers <- purrr::map(.x = kn_markers, .f = function(x) x$Marker)
kn_markers <- kn_markers[1:(length(kn_markers)-1)]

library(GSVA)
kn_markers_scores <- gsva(avg_exp, kn_markers)

library(ComplexHeatmap)

pdf('./figures/Cluster_GSVA_041421.pdf')
Heatmap(mtx_percentile(kn_markers_scores), name = 'Score')
dev.off()

# split group

Idents(test) <- test$group_cluster
avg_exp_group_cluster <- AverageExpression(test, assays = 'SCT')[[1]]

avg_exp_group_cluster <- avg_exp_group_cluster[,str_sort(colnames(avg_exp_group_cluster), numeric = T)]

library(circlize)



pdf('./figures/Cluster_GSVA_041421.pdf')
Heatmap(mtx_percentile(kn_markers_scores),
        col = colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027")),
        name = 'Score')
Heatmap(mtx_percentile(gsva(avg_exp_group_cluster[,str_detect(colnames(avg_exp_group_cluster), 'E')], 
                            kn_markers)),
        col = colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027")),
        name = 'Score',
        cluster_columns = F, column_title = 'ER') +
Heatmap(mtx_percentile(gsva(avg_exp_group_cluster[,str_detect(colnames(avg_exp_group_cluster), 'P')], 
                            kn_markers)),
        col = colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027")),
        name = 'Score',
        cluster_columns = F, column_title = 'PR')
dev.off()


# Go analysis

library(clusterProfiler)
library(org.Hs.eg.db)

cluster_go <- function(markers) {
    
    res <- list()
    
    n <- length(markers)
    
    for (i in 1:n) {
        
        df <- markers[[i]]
        go_res <- df %>%
            filter(p_val_adj <= 0.01 & avg_log2FC > 0.5)  %>%
            pull(Gene) %>%
            enrichGO(OrgDb = org.Hs.eg.db,
                     ont = 'BP',
                     keyType = 'SYMBOL',
                     qvalueCutoff = 0.1)
        
        go_res_df <- go_res@result %>%
            filter(p.adjust <= 0.05)
        
        res[[i]] <- go_res_df
    }
    return(res)
}

test_df <- markers[[1]]
test_df <- test_df %>%
    filter(p_val_adj <= 0.01 & avg_log2FC > 0.5)  %>%
    pull(Gene) %>%
    enrichGO(OrgDb = org.Hs.eg.db,
             ont = 'BP',
             keyType = 'SYMBOL',
             qvalueCutoff = 0.1)
View(test_df@result)

cluster_go_res <- cluster_go(markers)

cluster_go_res_plot <- map2(cluster_go_res, 
                           paste0('Cluster_',0:8),
                           function(x, y) {
    x %>%
        head(15) %>%
        ggplot(aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
        geom_bar(stat = 'identity') +
        coord_flip() +
        ggtitle(y) +
        theme(axis.title.y = element_blank())
})

pdf('./figures/Cluster_GO_031821.pdf', width = 15, height = 15)
wrap_plots(cluster_go_res_plot, ncol = 2)
dev.off()

