
library(tidyverse)

# Processing

md <- read_csv('./data/bulk_RNAseq/bulk_meta.csv', col_types = 'ccc')
md <- as.data.frame(md)
rownames(md) <- md$Sample
md$Response <- str_extract(md$Sample3, 'ER|PR')


rna_bulk <- list()

for (i in 1:length(list.files('./data/bulk_RNAseq/bulk'))) {
    rna_bulk[[i]] <- read_tsv(str_sort(list.files('./data/bulk_RNAseq/bulk', 
                                         full.names = T),
                              numeric = T)[i])
}

for (i in 1:length(rna_bulk)) {
    rna_bulk[[i]] <- rna_bulk[[i]] %>%
        select(gene_name, count) %>%
        drop_na() %>%
        unique() %>%
        mutate(sample = md$Sample[i])
}

rna_bulk <- do.call('rbind', rna_bulk)

rna_bulk <- rna_bulk %>%
    pivot_wider(names_from = sample, values_from = count, values_fn = median)

rna_bulk <- as.data.frame(rna_bulk)

rownames(rna_bulk) <- rna_bulk$gene_name

rna_bulk$gene_name <- NULL

rna_bulk <- as.matrix(rna_bulk)
rna_bulk <- round(rna_bulk)

# Analysis

library(DESeq2)

dds_bulk <- DESeqDataSetFromMatrix(countData = rna_bulk,
                                   colData = md,
                                   design= ~ Response)
dds_bulk <- DESeq(dds_bulk)
res_bulk <- results(dds_bulk, 
                    contrast = c("Response","PR","ER"))

plotMA(res_bulk, ylim = c(-6,6))

plotPCA(vst(dds_bulk, blind = FALSE), intgroup = 'Response')



