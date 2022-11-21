# Convert Bulk RNA Seq counts to TPM
#
# Abderrahman Day
#-------------------------------------------------------------------------------

library(tidyverse)
library(DGEobj.utils)
library(janitor)
library(biomaRt)

counts2tpm <- function(subread_counts_path, tpm_path){

# Get ensembl gene database to convert gene ensembl ID to gene symbol
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_symbols<-getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
                    mart = ensembl)

# Read data from subread_counts.txt (featurecounts output)
subread_counts <- read_table(subread_counts_path, skip = 1) 

# Extract gene counts matrix
counts <- subread_counts %>% 
  dplyr::select(Geneid, t1=PGDX8570T1_RNA, t2=PGDX8570T2_RNA) %>% 
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
  rownames_to_column(var = "geneid") %>% 
  left_join(gene_symbols, by = c("geneid" = "ensembl_gene_id")) %>% 
  dplyr::select(ensembl_gene_id = geneid,
                gene_symbol = hgnc_symbol,
                t1_tpm = t1,
                t2_tpm = t2)
tpm1$gene_symbol[tpm1$gene_symbol == ""] <- NA

# Write table
write_tsv(tpm1, file = tpm_path)

}

# Plots
#library(ggpubr)
#cpm_join %>% 
#  ggplot(aes(PGDX8570T1_RNA, logt1)) +
#  geom_point() +
#  geom_abline(slope = 1, intercept = 0) +
#  stat_cor() 
