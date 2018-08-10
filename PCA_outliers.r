library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")
library(biomaRt)
source("fun.R")

##################################
## Principal Component Analyses ##
##################################

# minus outliers, see which genes contribute the most to PC1

## Load gene data from CSV files ##
## ----------------------------- ##

CFZ_genes <- read.csv("CFZ_genes.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("coldata.csv", row.names = 1)

# From PCA plots and data from 'returnData' pca option
# outlier 6hr: CFZ-22-t6-r3
# outlier 24hr: CFZ-DMSO-t24-r2, CFZ-18-t24-r3, ?CFZ-26-t24-r1?
outlier_row_column_list <- c(which(rownames(coldata) == "CFZ-18-t24-r3"),which(rownames(coldata) == "CFZ-DMSO-t24-r2"),which(rownames(coldata) == "CFZ-22-t6-r3"))
CFZ_genes <- CFZ_genes[, -outlier_row_column_list]
coldata <- coldata[-outlier_row_column_list,]

# individual time point, subset of above data

CFZ_6hr <- CFZ_genes[,1:14]
coldata_6 <- coldata[1:14,]
CFZ_24hr <- CFZ_genes[,15:27]
coldata_24 <- coldata[15:27,]

## Using DESeq2 Package ##
## -------------------- ##

# 6 hour #
dds <-  DESeqDataSetFromMatrix(countData = CFZ_6hr,
                               colData = coldata_6,
                               design = ~ Compound)
dds

# Prefiltering, removing rows with sum of counts under 10 #
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Set the reference (control) level/ factor, against which others will be compared
dds$Compound <- relevel(dds$Compound, ref = "DMSO")

dds_deseq <- DESeq(dds)

# 24 hour # 
dds_24 <-  DESeqDataSetFromMatrix(countData = CFZ_24hr,
                                  colData = coldata_24,
                                  design = ~ Compound)
dds_24

## Pre-filtering ##
keep <- rowSums(counts(dds_24)) >= 10
dds_24 <- dds_24[keep,]
dds_24

# Set the reference (control) level/ factor, against which others will be compared
dds_24$Compound <- relevel(dds_24$Compound, ref = "DMSO")

# Contrasting compounds against DMSO, using p adjusted value of 5% as cut off.
dds_deseq24 <- DESeq(dds_24)

## Regularising transformation ##
## --------------------------- ##
# 6hr
rld <- rlog(dds_deseq, blind=FALSE)
# 24 hr 
rld24 <- rlog(dds_deseq24, blind=FALSE) 


## Clustering ##
## ---------- ##

## PCA Plot ##
plotPCA(rld, intgroup="Compound")
plotPCA(rld24, intgroup="Compound")

# PCA data #
pca_6 <- plotPCA(rld, intgroup="Compound", returnData = TRUE)
pca_24 <- plotPCA(rld24, intgroup="Compound", returnData = TRUE)

# Use function manually to see what genes contribute most to PC1 variation (pca rotation function)
# From both plots PC1 varaition seems to be between active and inactive/ control compounds. Therefore these are the genes that are most likely to be differentially expressed
rv <- rowVars(assay(rld))
rv24 <- rowVars(assay(rld24))

select_rv <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv)))]

select_rv24 <- order(rv24, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv24)))]
pca <- prcomp(t(assay(rld)[select_rv, ]))
pca24 <- prcomp(t(assay(rld24)[select_rv24, ]))

pca_rot <- as.data.frame(pca$rotation[,1:2]) # See PC1 and PC2
PC1_top <- pca_rot[rev(order(pca_rot$PC1)),] # Order by genes contributing most variation to PC1

pca_rot24 <- as.data.frame(pca24$rotation[,1:2]) # See PC1 and PC2
PC1_top24 <- pca_rot[rev(order(pca_rot24$PC1)),] # Order by genes contributing most variation to PC1

# Take the top 50 contributing genes from 6 hour and 24 hour and take the distinct genes
gene_join <- c()
for(i in 1:50){
  gene_join <- c(gene_join, rownames(PC1_top)[i], rownames(PC1_top24)[i])
}
string_genes <- unique(gene_join)[1:50]
#string_genes <- unique(c(rownames(PC1_top)[1:50], rownames(PC1_top24)[1:50])) # Bias towards 6 hour genes being nearer the top
genes_data <- gene_id_name_raw(string_genes)

## write csv for top 50 genes contributing to clustering ##
## ----------------------------------------------------- ##
write.csv(genes_data, file = "top50_PC1.csv", row.names = FALSE)

## Explore outliers ##
## ---------------- ##

# Look at variation along specific axis for outliers, see if genes contributing are involved in specific pathway, e.g. inflammatory
# Later...