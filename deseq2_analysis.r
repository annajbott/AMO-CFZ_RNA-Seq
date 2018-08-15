library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")
source("fun.R")
library(stringr)
library(ggplot2)

#####################
## DESeq2 Analysis ##
#####################

## Load gene data from CSV files ##
## ----------------------------- ##

CFZ_genes <- read.csv("CFZ_genes.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("coldata.csv", row.names = 1)

# individual time point loads
CFZ_6hr <- read.csv("CFZ_genes_6hr.csv", row.names = 1, check.names = FALSE)
coldata_6 <- read.csv("coldata_6hr.csv", row.names = 1)
CFZ_24hr <- read.csv("CFZ_genes_24hr.csv", row.names = 1, check.names = FALSE)
coldata_24 <- read.csv("coldata_24hr.csv", row.names = 1)

## Set-up dds data from DESeq function ##
## ----------------------------------- ##

# 6hr data first #
# -------------- #

dds <-  DESeqDataSetFromMatrix(countData = CFZ_6hr,
                               colData = coldata_6,
                               design = ~ Compound)
dds

## Prefiltering, conservatively removing rows with under 10 reads ##

# Keep only rows with reads over 10 counts, so that memory of object is reduced
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Set the reference (control) level/ factor, against which others will be compared
dds$Compound <- relevel(dds$Compound, ref = "DMSO")

## Collapse technical replicates??? (Ask adam) ##
# If you collapse technical replicates will we lose all our variation?

# Contrasting compounds against DMSO, using p adjusted value of 5% as cut off.
dds_deseq <- DESeq(dds)

## 24 hour ##
## ------- ##

dds_24 <-  DESeqDataSetFromMatrix(countData = CFZ_24hr,
                                  colData = coldata_24,
                                  design = ~ Compound)
dds_24

## Pre-filtering ##
# Keep only rows with reads over 10 counts, so that memory of object is reduced
keep <- rowSums(counts(dds_24)) >= 10
dds_24 <- dds_24[keep,]
dds_24

# Set the reference (control) level/ factor, against which others will be compared
dds_24$Compound <- relevel(dds_24$Compound, ref = "DMSO")

# Contrasting compounds against DMSO, using p adjusted value of 5% as cut off.
dds_deseq24 <- DESeq(dds_24)

######################################
## Differential Expression Analysis ##
######################################

# 6 hour results, order padj all under 0.05 in order of magnitude (descending) of log2 fold change.
res_26_lfc <- results_process(dds_deseq, contrast = c("Compound", "26", "DMSO"), alpha = 0.05 )
res_22_lfc <- results_process(dds_deseq, contrast = c("Compound", "22", "DMSO"), alpha = 0.05 )
res_13_lfc <- results_process(dds_deseq, contrast = c("Compound", "13", "DMSO"), alpha = 0.05 )
res_18_lfc <- results_process(dds_deseq, contrast = c("Compound", "18", "DMSO"), alpha = 0.05 )
# NCP 22 has no genes with adj p value under 0.05.

## 24 hour results ##
res_26_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "26", "DMSO"), alpha = 0.05 )
res_22_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "22", "DMSO"), alpha = 0.05 )
res_13_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "13", "DMSO"), alpha = 0.05 )
res_18_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "18", "DMSO"), alpha = 0.05 )
# NCP 22 and MAZ 18 compounds have no genes with a padj value <0.05, compared to DMSO. 

# All unique differentially expressed genes (padj <0.05) with accompanying gene hgnc symbols, for 6hr and 24 hr
#26
genes_dif_expressed_26_6 <- gene_id_name(res_26_lfc, id_name_table = TRUE, top = nrow(res_26_lfc))
genes_dif_expressed_26_24 <- gene_id_name(res_26_lfc_24, id_name_table = TRUE, top = nrow(res_26_lfc_24))
genes_dif_expressed_26 <- unique(rbind(genes_dif_expressed_26_6[,c(1,2,4)],genes_dif_expressed_26_24[,c(1,2,4)]))
# 13
genes_dif_expressed_13_6 <- gene_id_name(res_13_lfc, id_name_table = TRUE, top = nrow(res_13_lfc))
genes_dif_expressed_13_24 <- gene_id_name(res_13_lfc_24, id_name_table = TRUE, top = nrow(res_13_lfc_24))
genes_dif_expressed_13 <- unique(rbind(genes_dif_expressed_13_6[,c(1,2,4)],genes_dif_expressed_13_24[,c(1,2,4)]))


# Top 50 genes (greatest magnitude fold change)
top50_res_26_lfc <- gene_id_name(res_26_lfc, id_name_table = TRUE)
top50_res_26_lfc$Compound <- "26"
top50_res_26_lfc$TimePoint <- 6
top50_res_13_lfc <- gene_id_name(res_13_lfc, id_name_table = TRUE)
top50_res_13_lfc$Compound <- "13"
top50_res_13_lfc$TimePoint <- 6
top50_res_18_lfc <- gene_id_name(res_18_lfc, id_name_table = TRUE)
top50_res_18_lfc$Compound <- "18"
top50_res_18_lfc$TimePoint <- 6
# 24
top50_res_26_lfc_24 <- gene_id_name(res_26_lfc_24, id_name_table = TRUE)
top50_res_26_lfc_24$Compound <- "26"
top50_res_26_lfc_24$TimePoint <- 24
top50_res_13_lfc_24 <- gene_id_name(res_13_lfc_24, id_name_table = TRUE)
top50_res_13_lfc_24$Compound <- "13"
top50_res_13_lfc_24$TimePoint <- 24

top50_all <- rbind(top50_res_26_lfc,top50_res_13_lfc, top50_res_18_lfc,top50_res_26_lfc_24,top50_res_13_lfc_24)
top50_all <- top50_all[!(is.na(top50_all$hgnc_symbol) | top50_all$hgnc_symbol==""), ]

# Write csv files for genes #
# ------------------------- #

write.csv(top50_all, file = "top50lfc_difexpressed_all_samples.csv", row.names = FALSE)

#############
### Plots ###
#############

# Assorted mess #
## MA-plots ##
# 6 hour
res_26_noshrink <- results(dds_deseq, contrast = c("Compound", "26","DMSO"), alpha = 0.05)
res_22_noshrink <- results(dds_deseq, contrast = c("Compound", "22","DMSO"), alpha = 0.05)
res_13_noshrink <- results(dds_deseq, contrast = c("Compound", "13","DMSO"), alpha = 0.05)
res_18_noshrink <- results(dds_deseq, contrast = c("Compound", "18","DMSO"), alpha = 0.05)

plotMA(res_26_noshrink, ylim = c(-2,2))
plotMA(res_22_noshrink, ylim = c(-2,2))
plotMA(res_13_noshrink, ylim = c(-2,2))
plotMA(res_18_noshrink, ylim = c(-2,2))
# 24 hour
plotMA(res_26_lfc_24, ylim = c(-2,2))
plotMA(res_22_lfc_24, ylim = c(-2,2))
plotMA(res_13_lfc_24, ylim = c(-2,2))
plotMA(res_18_lfc_24, ylim = c(-2,2))

# Apply a regularised log transformation. Variance stabilising effect. Similar to VST
# Useful in checking for outliers. Useful for visualisation, clistering or ML.

rld <- rlog(dds_deseq, blind=FALSE)
# 24
rld24 <- rlog(dds_deseq24, blind=FALSE) 
# Blind true, compare samples unbiased by prior info on samples. False for downstream analysis, maybe should be true?

## PCA plots ##
plotPCA(rld, intgroup="Compound")
plotPCA(rld24, intgroup="Compound")

# 2 outliers from 24 hour PCA, 1 DMSO and 1 MAZ 18. Affects PC1
rv <- rowVars(assay(rld24))
select_rv <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv)))]
pca <- prcomp(t(assay(rld24)[select_rv, ]))

pca_rot <- as.data.frame(pca$rotation[,1:2])
PC1 <- pca_rot[rev(order(pca_rot$PC1)),][1:30,] # Order by genes contributing most variation to PC1
PC1 <- cbind(gene_id_name_raw(rownames(PC1)),PC1)[,1:3]


# Heat map, sample distances #

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Heat map 24
sampleDists <- dist(t(assay(rld24)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Plots of AAR known genes ##
## ------------------------ ##

samples_list <- list(genes_dif_expressed_13_6, genes_dif_expressed_13_24, genes_dif_expressed_26_6, genes_dif_expressed_26_24)
gene_list <- c("EIF4EBP1", "ASNS", "GPT2", "DDIT3" , "ATF4", "SLC6A9")
plot_table <- data.frame()

plot_table[1:2,1] <- "13"
plot_table[3:4,1] <- "26"
plot_table[c(1,3),2] <- "6"
plot_table[c(2,4),2] <- "24"
rownames(plot_table) <- c("13_6", "13_24", "26_6", "26_24")
colnames(plot_table) <- c("Compound", "TimePoint")

for(i in seq_along(samples_list)){
  for(j in seq_along(gene_list)){
    subset_data <- samples_list[[i]] %>% filter(str_detect(hgnc_symbol, gene_list[j]))
    plot_table[i,j + 2] <- subset_data[1,"log2FoldChange"]
  }
}
colnames(plot_table) <- c("Compound", "TimePoint",gene_list)
plot_table$Compound <- as.factor(plot_table$Compound)
plot_table$TimePoint <- as.factor(plot_table$TimePoint)

plot_table_tidy <- gather(plot_table,"gene", "Log2FoldChange", -c(Compound, TimePoint))

ggplot(data=plot_table_tidy, aes(x=gene, y=Log2FoldChange, fill = Compound)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge())

## 6hr vs 24 hr plots ##
## ------------------ ##

## Plot of gene log fold change 6 hr vs 24 hr ##
all_results_26_6 <- results_process(dds_deseq, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE)
all_results_26_24 <- results_process(dds_deseq24, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE)
all_results_13_6 <- results_process(dds_deseq, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE)
all_results_13_24 <- results_process(dds_deseq24, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE)

# Get rid of NAs and blank p values. Make match
#26
all_results_26_6 <- all_results_26_6[!(is.na(all_results_26_6$padj) | all_results_26_6$padj==""), ]
all_results_26_24 <- all_results_26_24[!(is.na(all_results_26_24$padj) | all_results_26_24$padj==""), ]
#13
all_results_13_6 <- all_results_13_6[!(is.na(all_results_13_6$padj) | all_results_13_6$padj==""), ]
all_results_13_24 <- all_results_13_24[!(is.na(all_results_13_24$padj) | all_results_13_24$padj==""), ]

# 26 keep
kept <- all_results_26_6[rownames(all_results_26_6) %in% rownames(all_results_26_24),]
all(rownames(all_results_26_24) %in% rownames(kept))
kept2 <- all_results_26_24[rownames(all_results_26_24) %in% rownames(kept),]
# check in the right order
all(rownames(kept) == rownames(kept2)) # Right order

# 13 keep
kept_13 <- all_results_13_6[rownames(all_results_13_6) %in% rownames(all_results_13_24),]
all(rownames(all_results_13_24) %in% rownames(kept))
kept2_13 <- all_results_13_24[rownames(all_results_13_24) %in% rownames(kept_13),]
# check in the right order
all(rownames(kept_13) == rownames(kept2_13)) # Right order

# 26 data frame and plot
time_vs_26 <- as.data.frame(cbind(kept$log2FoldChange, kept2$log2FoldChange))
rownames(time_vs_26) <- rownames(kept)
colnames(time_vs_26) <- c("log_fold_change_6hr", "log_fold_change_24hr")
time_vs_26$threshold_status <- ifelse(time_vs_26$log_fold_change_6hr > 1 & time_vs_26$log_fold_change_24hr > 1, TRUE, FALSE )
ggplot(time_vs_26, aes(x=log_fold_change_6hr, y=log_fold_change_24hr, color = threshold_status)) + geom_point() +
 labs(x = "Log 2 Fold Change (NCP 26 vs DMSO)- 6hr", y =  "Log 2 Fold Change (NCP 26 vs DMSO)- 24hr" ) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

# 13 data frame and plot
time_vs_13 <- as.data.frame(cbind(kept_13$log2FoldChange, kept2_13$log2FoldChange))
rownames(time_vs_13) <- rownames(kept_13)
colnames(time_vs_13) <- c("log_fold_change_6hr", "log_fold_change_24hr")
# New column indicating if gene expression has doubled for both 6hr and 24 hr, can then colour graph using this variable
time_vs_13$threshold_status <- ifelse(time_vs_13$log_fold_change_6hr > 1 & time_vs_13$log_fold_change_24hr > 1, TRUE, FALSE )
ggplot(time_vs_13, aes(x=log_fold_change_6hr, y=log_fold_change_24hr, color = threshold_status)) + 
  geom_point() + 
  labs(x = "Log 2 Fold Change (MAZ 13 vs DMSO)- 6hr", y =  "Log 2 Fold Change (MAZ 13 vs DMSO)- 24hr" ) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
