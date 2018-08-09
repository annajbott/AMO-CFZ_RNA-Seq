library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")
library(biomaRt)
source("fun.R")

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

# Add gene names to tables
top50_res_26_lfc <- gene_id_name(res_26_lfc)
top50_res_13_lfc <- gene_id_name(res_13_lfc)
top50_res_18_lfc <- gene_id_name(res_18_lfc)
# 24
top50_res_26_lfc_24 <- gene_id_name(res_26_lfc_24)
top50_res_13_lfc_24 <- gene_id_name(res_13_lfc_24)

##

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

# test for getting PCA data
rv <- rowVars(assay(rld))
select_rv <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv)))]
pca <- prcomp(t(assay(rld)[select_rv, ]))

pca_rot <- as.data.frame(pca$rotation[,1:2])
pca_PC2 <- pca_rot[rev(order(pca_rot$PC2)),] # Order by genes contributing most variation to PC2

pca_6 <- plotPCA(rld, intgroup=c("Compound"), returnData = TRUE)
# CFZ-22-t6-r3 seems to be an outlier

######## Investigate PCA outlier ######## 

genes_22 <- select(CFZ_6hr, contains("-22-"))
keep <- rowSums(genes_22) >= 20
genes_22 <- genes_22[keep,]
genes_22 <- filter(genes_22, !is.na(genes_22))

pca_PC2_large <- pca_PC2[1:50,]
genes_22_pc2 <- genes_22[row.names(pca_PC2_large),]

genes_22_pc2 <- transform(genes_22_pc2, percent= abs(`CFZ-22-t6-r3` - mean(c(`CFZ-22-t6-r1`, `CFZ-22-t6-r2`)))/mean(c(`CFZ-22-t6-r1`, `CFZ-22-t6-r2`)), check.names = FALSE)
genes_22_pc2

########
# Heat map, explore clusters
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



rv <- rowVars(assay(rld24))
select_rv <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv)))]
pca <- prcomp(t(assay(rld24)[select_rv, ]))

pca_rot <- as.data.frame(pca$rotation[,1:2])
#pca_PC2 <- pca_rot[rev(order(pca_rot$PC2)),] # Order by genes contributing most variation to PC2

pca_24 <- plotPCA(rld24, intgroup=c("Compound"), returnData = TRUE)
