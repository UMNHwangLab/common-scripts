# Convert Bulk RNA Seq counts to TPM
#
# Abderrahman Day
#-------------------------------------------------------------------------------

library(tidyverse)
library(DGEobj.utils)
library(janitor)
library(biomaRt)

counts2tpm <- function(counts_path, species = "human"){
  
  # Get ensemble and gene symbol databases

  if (species == "human"){
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    attributes = c('hgnc_symbol','ensembl_gene_id')
    gene_symbols<-getBM(attributes=attributes, mart = ensembl) %>% 
      dplyr::rename(ensembl_gene_id = 2, gene_symbol = 1)
  }
  else if (species == "mouse"){
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    attributes <- c("ensembl_gene_id", "mgi_symbol")
    gene_symbols<-getBM(attributes=attributes, mart = ensembl) %>% 
      dplyr::rename(ensembl_gene_id = 1, gene_symbol = 2)
  }
  else{
  }
  
  # Get ensembl and gene symbol for MOUSE
  
  # Read data from subread_counts.txt (featurecounts output)
  subread_counts <- read_table(counts_path, skip = 1) 
  
  # Extract gene counts matrix
  counts <- subread_counts %>% 
    dplyr::select(-c(2:6)) %>% 
    column_to_rownames(var = "Geneid") %>% 
    as.matrix()
  
  # Extract gene lengths matrix
  lengths <- subread_counts %>% 
    dplyr::select(Geneid, length = Length) %>% 
    column_to_rownames(var = "Geneid") %>% 
    as.matrix()
  
  # Calculate CPM
  #cpm1 <- convertCounts(counts, "CPM") %>% 
  #  as.data.frame() %>% 
  #  rownames_to_column(var = "geneid")
  
  # Calculate TPM, add gene symbols, convert empty gene symbols to NA
  tpm1 <- convertCounts(counts, "TPM", geneLength = lengths) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ensembl_gene_id") %>% 
    left_join(gene_symbols, by = "ensembl_gene_id") %>% 
    relocate(ensembl_gene_id,gene_symbol)
  
  tpm1$gene_symbol[tpm1$gene_symbol == ""] <- NA

 return(tpm1) 
}

# Plots
#library(ggpubr)
#cpm_join %>% 
#  ggplot(aes(PGDX8570T1_RNA, logt1)) +
#  geom_point() +
#  geom_abline(slope = 1, intercept = 0) +
#  stat_cor() 
