
library(tidyverse)

df_lr <- read_csv("./data/ligand_receptor_list.csv")

lr_score <- function(mtx, cluster_1, cluster_2, ligand, receptor) {
    
    if (!ligand %in% rownames(mtx) | !receptor %in% rownames(mtx)) {
        return(0)
    } else {
        m <- mean(mtx)
        p <- mtx[ligand, cluster_1] * mtx[receptor, cluster_2]
        score <- sqrt(p) / (m + sqrt(p))
        return(score)
    }
}

df_lr_score <- expand_grid(df_lr, colnames(avg_exp))
colnames(df_lr_score)[3] <- 'Cluster_1'

df_lr_score <- expand_grid(df_lr_score, colnames(avg_exp))
colnames(df_lr_score)[4] <- 'Cluster_2'


df_lr_score$Score <- purrr::pmap_dbl(.l = list(df_lr_score$Cluster_1,
                                               df_lr_score$Cluster_2,
                                               df_lr_score$Ligand,
                                               df_lr_score$Receptor), .f = lr_score, mtx = avg_exp)

df_lr_score$Pair <- paste0(df_lr_score$Ligand, '_', df_lr_score$Receptor)
df_lr_score$Cluster <- paste0(df_lr_score$Cluster_1, '_', df_lr_score$Cluster_2)


mtx_lr_score <- df_lr_score %>%
    select(Score:Cluster) %>%
    pivot_wider(names_from = 'Cluster', values_from = "Score") %>%
    as.data.frame()

rownames(mtx_lr_score) <- mtx_lr_score$Pair
mtx_lr_score$Pair <- NULL
mtx_lr_score <- as.matrix(mtx_lr_score)

library(ComplexHeatmap)
library(circlize)

pdf('./figures/Cluster_LR_053021.pdf', width = 8.5, height = 11)
Heatmap(mtx_lr_score[rowSums(mtx_lr_score) > 60, 
                     c("Cluster_5_Cluster_2", "Cluster_2_Cluster_5",
                       "Cluster_5_Cluster_1", "Cluster_1_Cluster_5",
                       "Cluster_3_Cluster_8", "Cluster_8_Cluster_3",
                       "Cluster_8_Cluster_2", "Cluster_2_Cluster_8",
                       "Cluster_8_Cluster_4", "Cluster_4_Cluster_8",
                       "Cluster_3_Cluster_4", "Cluster_4_Cluster_3",
                       "Cluster_4_Cluster_0", "Cluster_0_Cluster_4")], name = 'Score',
        column_split = rep(1:7, each = 2),
        col = colorRamp2(c(0.65, 0.85, 0.95), c("#4575b4", "white", "#d73027")), cluster_columns = F)
dev.off()

mtx_lr_score[rowSums(mtx_lr_score) > 60,] %>% dim()

Heatmap(mtx_lr_score[rowSums(mtx_lr_score) > 60,], name = 'Score',
        col = colorRamp2(c(0.5, 0.85, 0.95), c("#4575b4", "white", "#d73027")), cluster_columns = F)

Heatmap(mtx_lr_score)

rep(c("A", "B"), each = 9)

# 060221

# ER

df_lr_score_er <- expand_grid(df_lr, colnames(ov_avg_exp_group_cluster[,1:8]))
colnames(df_lr_score_er)[3] <- 'Cluster_1'

df_lr_score_er <- expand_grid(df_lr_score_er, colnames(ov_avg_exp_group_cluster[,1:8]))
colnames(df_lr_score_er)[4] <- 'Cluster_2'

df_lr_score_er$Score <- purrr::pmap_dbl(.l = list(df_lr_score_er$Cluster_1,
                                                  df_lr_score_er$Cluster_2,
                                                  df_lr_score_er$Ligand,
                                                  df_lr_score_er$Receptor), 
                                        .f = lr_score, 
                                        mtx = ov_avg_exp_group_cluster[,1:8])

df_lr_score_er$Pair <- paste0(df_lr_score_er$Ligand, '_', df_lr_score_er$Receptor)
df_lr_score_er$Cluster <- paste0(df_lr_score_er$Cluster_1, '_', df_lr_score_er$Cluster_2)

mtx_lr_score_er <- df_lr_score_er %>%
    select(Score:Cluster) %>%
    pivot_wider(names_from = 'Cluster', values_from = "Score") %>%
    as.data.frame()

rownames(mtx_lr_score_er) <- mtx_lr_score_er$Pair
mtx_lr_score_er$Pair <- NULL
mtx_lr_score_er <- as.matrix(mtx_lr_score_er)

# PR

df_lr_score_pr <- expand_grid(df_lr, colnames(ov_avg_exp_group_cluster[,9:16]))
colnames(df_lr_score_pr)[3] <- 'Cluster_1'

df_lr_score_pr <- expand_grid(df_lr_score_pr, colnames(ov_avg_exp_group_cluster[,9:16]))
colnames(df_lr_score_pr)[4] <- 'Cluster_2'

df_lr_score_pr$Score <- purrr::pmap_dbl(.l = list(df_lr_score_pr$Cluster_1,
                                                  df_lr_score_pr$Cluster_2,
                                                  df_lr_score_pr$Ligand,
                                                  df_lr_score_pr$Receptor), 
                                        .f = lr_score, 
                                        mtx = ov_avg_exp_group_cluster[,9:16])

df_lr_score_pr$Pair <- paste0(df_lr_score_pr$Ligand, '_', df_lr_score_pr$Receptor)
df_lr_score_pr$Cluster <- paste0(df_lr_score_pr$Cluster_1, '_', df_lr_score_pr$Cluster_2)

mtx_lr_score_pr <- df_lr_score_pr %>%
    select(Score:Cluster) %>%
    pivot_wider(names_from = 'Cluster', values_from = "Score") %>%
    as.data.frame()

rownames(mtx_lr_score_pr) <- mtx_lr_score_pr$Pair
mtx_lr_score_pr$Pair <- NULL
mtx_lr_score_pr <- as.matrix(mtx_lr_score_pr)


pdf('./figures/Cluster_LR_060221.pdf', width = 8.5, height = 11)
Heatmap(mtx_lr_score[rowSums(mtx_lr_score) > 60, 
                     c("Cluster_5_Cluster_2", "Cluster_2_Cluster_5",
                       "Cluster_8_Cluster_1", "Cluster_1_Cluster_8",
                       "Cluster_8_Cluster_2", "Cluster_2_Cluster_8",
                       "Cluster_8_Cluster_3", "Cluster_3_Cluster_8",
                       "Cluster_8_Cluster_4", "Cluster_4_Cluster_8")], name = 'Score',
        column_split = rep(1:5, each = 2),
        col = colorRamp2(c(0.65, 0.8, 0.95), c("#4575b4", "white", "#d73027")), cluster_columns = F, 
        column_title = 'ALL')
Heatmap(mtx_lr_score_er[rowSums(mtx_lr_score_er) > 50, 
                     c("ER_5_ER_2", "ER_2_ER_5",
                       "ER_8_ER_1", "ER_1_ER_8",
                       "ER_8_ER_2", "ER_2_ER_8",
                       "ER_8_ER_3", "ER_3_ER_8",
                       "ER_8_ER_4", "ER_4_ER_8")], name = 'Score',
        column_split = rep(1:5, each = 2),
        col = colorRamp2(c(0.7, 0.85, 0.95), c("#4575b4", "white", "#d73027")), cluster_columns = F,
        column_title = 'ER')
Heatmap(mtx_lr_score_pr[rowSums(mtx_lr_score_pr) > 45, 
                        c("PR_5_PR_2", "PR_2_PR_5",
                          "PR_8_PR_1", "PR_1_PR_8",
                          "PR_8_PR_2", "PR_2_PR_8",
                          "PR_8_PR_3", "PR_3_PR_8",
                          "PR_8_PR_4", "PR_4_PR_8")], name = 'Score',
        column_split = rep(1:5, each = 2),
        col = colorRamp2(c(0, 0.75, 0.95), c("#4575b4", "white", "#d73027")), cluster_columns = F,
        column_title = 'PR')
dev.off()


# output score

library(openxlsx)

openxlsx::write.xlsx(list(mtx_lr_score, mtx_lr_score_er, mtx_lr_score_pr),
                     "./tables/LR_Scores_060221.xlsx", row.names = TRUE)




