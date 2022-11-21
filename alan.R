# ALAN Function
#
# Purpose: This script takes gene expression data from multiple samples and returns 
#          the associations between genes 
#
# Required inputs: Dataframe of gene expression data
#                   -Each row represents a gene
#                   -Each column represents a sample
#                   -First column contains the gene names
#                   -First row contains sample names
#
# Output: (1) TSV file of ALAN export matrix 1 (gene correlations)
#         (2) TSV file of ALAN export matrix 2 (gene network correlations)
#	  (3) Returns Dataframe of ALAB export matrix 2
#
# Abderrahman Day, Ashraf Shabaneh, and Hannah Bergom
# Last updated 2022-10-27 by A.D. 
#-------------------------------------------------------------------------------

alan_fun <- function(df, alan1_path, alan2_path, delim_out){
  
  
  alan1 <- df %>%                                           # Export matrix 1       
    column_to_rownames(var=colnames(.)[1]) %>%                # convert 1st column to rownames
    na.omit() %>%                                             # omit rows (genes) with any NAs
    t() %>%                                                   # transpose so that 1 column = 1 gene
    as.data.frame() %>%                                       # convert to df from matrix
    .[,sapply(., function(v) var(v, na.rm=TRUE)!=0)] %>%      # remove columns (genes) with zero variance
    cor(method = "spearman")                                  # Spearman correlation

  alan1 %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    write_delim(alan1_path, delim = delim_out)
  
  alan2 <- alan1 %>% 
    cor(method = "pearson")                                   # Pearson correlation 

  alan2  %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>%
    write_delim(alan2_path, delim = delim_out)

  return(alan2)

}
