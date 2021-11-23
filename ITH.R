
library(DEPTH)
library(tidyverse)


DEPTH <- function(exp, match) {
    #Two files need to be input into this function.
    exp <- as.matrix(exp)
    samp <- colnames(exp)
    gene <- rownames(exp)
    
    match <- as.matrix(match)
    #Pick up normal samples.
    nor_pos <- which(samp %in% match[which(match[, 2] == "Normal"), 1]) 
    #Pick up tumor samples.
    tum_pos <- which(samp %in% match[which(match[, 2] == "Tumor"), 1]) 
    exp_tum <- exp[, tum_pos]; samp_tum <- samp[tum_pos]; 
    score <- matrix(0, nrow <- dim(exp_tum)[1], ncol <- dim(exp_tum)[2])
    
    if(length(nor_pos) > 0){
        nor <- c(); 
        for(j in 1:dim(exp)[1]){
            nor[j] <- mean(as.numeric(exp[j, nor_pos]))
        }
        #Calculate the average values of each gene in normal sample.
        #Calculate the heterogeneity score of each gene.
        for(s in 1 : dim(exp_tum)[1]){
            for(u in 1 : dim(exp_tum)[2]){
                score[s, u] <- (as.numeric(exp_tum[s, u]) - as.numeric(nor[s]))^2
            }
        }
    }else if(length(nor_pos)==0){
        for(s in 1 : dim(exp_tum)[1]){ 
            for(u in 1 : dim(exp_tum)[2]){
                score[s, u] <- (as.numeric(exp_tum[s, u]) - mean(as.numeric(exp_tum[s, ])))^2
            }
        }
    }
    colnames(score) <- samp_tum; rownames(score) <- gene;
    heterogeneity_score <- c();
    for(z in 1:length(samp_tum)){
        heterogeneity_score[z] <- sd(as.numeric(score[, z]))
    }
    heterogeneity_score <- cbind(samp_tum, heterogeneity_score); 
    #calculate the heterogeneity score of each sample.
    colnames(heterogeneity_score) <- c("sample", "score")
    #DEPTH function will output the heterogeneity score of each tumor sample.
    heterogeneity_score <- as.data.frame(heterogeneity_score)
    heterogeneity_score$score <- as.numeric(heterogeneity_score$score)
    rownames(heterogeneity_score) <- NULL
    return(heterogeneity_score)
}


# calculation

ith_cluster <- DEPTH(ov_avg_exp_cluster,
                     tibble(Sample = colnames(ov_avg_exp_cluster),
                            Identification = 'Tumor'))

ith_group_cluster <- DEPTH(ov_avg_exp_group_cluster,
                           tibble(Sample = colnames(ov_avg_exp_group_cluster),
                                  Identification = 'Tumor'))

ith_bulk <- DEPTH(log2(rna_bulk + 1),
                  tibble(Sample = colnames(rna_bulk),
                         Identification = 'Tumor'))

ith_group_cluster <- ith_group_cluster %>%
    mutate(group = str_extract(sample,'(ER|PR)'))



# plot
pdf('./figures/ITH_052621.pdf')

ith_bulk %>%
    mutate(group = rep(c('ER','PR'), each = 6)) %>%
    ggplot(aes(x = group, y = score, fill = group)) +
    geom_boxplot() +
    geom_point() +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw() +
    ggtitle('Bulk')

ith_group_cluster %>%
    ggplot() +
    geom_bar(aes(x = sample, y = score*10, fill = group), stat = 'identity') +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw() +
    ggtitle('Spatial') +
    scale_y_log10()
    

dev.off()





