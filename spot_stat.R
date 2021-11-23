
library(tidyverse)
library(Seurat)

# add group

test@meta.data$group <- ifelse(str_detect(test@meta.data$orig.ident, 'E'), 'ER', 'PR')


pdf('./figures/Cluster_Size_041421.pdf', width = 8.5, height = 7)

test@meta.data %>%
    as_tibble() %>%
    ggplot(aes(x = seurat_clusters, fill = group)) +
    geom_bar(position = 'dodge') +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw() +
    ggtitle('Cluster Size between ER and PR - No. of Spots')

test@meta.data %>%
    as_tibble() %>%
    ggplot(aes(x = seurat_clusters, fill = group)) +
    geom_bar(position = 'fill') +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw() +
    ggtitle('Cluster Size between ER and PR - Proportion')

dev.off()
