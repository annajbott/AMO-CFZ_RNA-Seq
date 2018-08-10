library(tidyverse)
source("fun.R")

#########################
## Gene cross checking ##
#########################

# Top 50 genes contributing the most to PC1 variation (difference between active and inactive compounds)
PC1_top_genes <- read_csv(file = "top50_PC1.csv")
# genes with largest log fold change, (p<0.05)
dif_expressed_top_all <- read_csv(file = "top50lfc_difexpressed_all_samples.csv")
dif_expressed_top_all$Compound <- as.factor(dif_expressed_top_all$Compound)

# Only keep 13 and 26 samples
dif_expressed_top_26_13 <- filter(dif_expressed_top_all, Compound != "18")
dif_expressed_top_26_13 <- dplyr::select(dif_expressed_top_26_13, -c(Compound, TimePoint))

# Keep distinct genes, no repeats
dif_expressed_top_26_13<- unique(dif_expressed_top_26_13)
# Cross check top sample genes with cluster PC1 variation contributing genes. Keep those in both
keep <- dif_expressed_top_26_13$ensembl_gene_id %in% PC1_top_genes$ensembl_gene_id
dif_expressed_top_26_13 <- dif_expressed_top_26_13[keep,]

# Cross check top sample genes with cluster PC1 variation contributing genes. Concatanate lists
dif_expressed_top_26_13$ensembl_gene_id %in% PC1_top_genes$ensembl_gene_id
genes_together <- unique(rbind(dif_expressed_top_26_13,PC1_top_genes))

which(PC1_top_genes$ensembl_gene_id %!in% dif_expressed_top_26_13$ensembl_gene_id)

genes_together
