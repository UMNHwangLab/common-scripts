# alan_geneList_score.R
#
# ALAN Gene List Scoring Function.
#
# Usage: alan_geneList_score(x, gene_list, gene_list_name = "gene_list")
#
# INPUTS: [df] (n x 2) - columns: [gene, ALAN_score]
#         [list] - Gene list
#         [character] - name of gene list
          
# OUTPUTS:
#         -Plot of ROC curve
#         -AUC
#         -Mean ALAN score of gene list
#         -Mean ALAN rank of gene list

# Load Packages
library(tidyverse)
library(pROC)

# ALAN GENE LIST SCORING FUNCTION
alan_geneList_score <- function(x, gene_list, gene_list_name = "gene_list"){
  # AR NELSON ALAN
  gene_name <- colnames(x[2])
  data <- x %>% 
    rename(id = 1, alan_score = 2) %>% 
    arrange(desc(alan_score)) %>% 
    mutate(in_geneList = (id %in% gene_list),
           cumSum = cumsum(in_geneList),
           cumFrac = 100*cumSum / length(gene_list),
           rank = row_number(),
           rank_fraction = 100*rank / nrow(.))
  
  p1 <- data %>% 
    ggplot(aes(rank_fraction, cumFrac)) +
    coord_fixed() +
    geom_line() +
    geom_abline(slope = 1, color = "red") +
    annotate("text", x = 49, y = 52, label = "Random Ranks Performance", angle='45', color = "red") +
    theme_bw() +
    scale_x_continuous(limits=c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0, 100), expand = c(0, 0)) +
    labs(title = sprintf("Ranking Accuracy of ALAN for %s gene list (%d genes)", gene_list_name, length(gene_list)),
         x = sprintf("%s ALAN Score Rank (percent of total)", gene_name),
         y = "Cumulative fraction of gene list captured (percent)")
  
  gene_list_ranks <- data.frame(id = gene_list) %>% 
    left_join(data, by = "id") 
  
  # Stats
  n_overlapping_genes = length(intersect(gene_list, data$id))
  pct_overlapping_genes = round(100*n_overlapping_genes/length(gene_list),1)
  mean_rank = mean(gene_list_ranks$rank, na.rm = TRUE)
  mean_rank_neat =  format(round(mean_rank,1), nsmall=1, big.mark=",")
  mean_rank_pct = round(100 * mean_rank / nrow(data),0)
  median_rank = median(gene_list_ranks$rank, na.rm = TRUE)
  median_rank_neat = format(median_rank, nsmall=1, big.mark=",")
  median_rank_pct = round(100 * median_rank / nrow(data),0)
  mean_alan_score = mean(gene_list_ranks$alan_score, na.rm = TRUE)
  roc_obj <- roc(data$in_geneList, data$rank)
  auc = roc_obj$auc[[1]]
  
  
  results <- paste0(
    sprintf("Ranking Accuracy of ALAN for %s gene list (%d genes)\n", gene_list_name, length(gene_list)),
    "-------------------------------------------------------------\n",
    sprintf("Number genes present in gene list and ALAN dataset: %d (%.1f%%)\n", n_overlapping_genes, pct_overlapping_genes),
    sprintf("Mean ALAN Rank of %d overlapping genes: %s (top %.0f%% of genes)\n", n_overlapping_genes, mean_rank_neat, mean_rank_pct),
    sprintf("Median ALAN Rank of %d overlapping genes: %s (top %.0f%% of genes)\n", n_overlapping_genes, median_rank_neat, median_rank_pct),
    sprintf("Mean ALAN Score of %d overlapping genes: %.3f\n", n_overlapping_genes, mean_alan_score),
    sprintf("Area under ROC curve: %.3f\n",auc)
  )
  
  writeLines(results)
  
  
  return(p1)
}