library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")

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
results_26 <- results(dds_deseq, contrast = c("Compound", "26", "DMSO"), alpha = 0.05)
results_22 <- results(dds_deseq, contrast = c("Compound", "22", "DMSO"), alpha = 0.05)
results_13 <- results(dds_deseq, contrast = c("Compound", "13", "DMSO"), alpha = 0.05)
results_18 <- results(dds_deseq, contrast = c("Compound", "18", "DMSO"), alpha = 0.05)
# Results with p values over 5%
results_26_p <- results_26[!is.na(results_26$padj) & results_26$padj<= 0.05,]
results_22_p <- results_22[!is.na(results_22$padj) & results_22$padj<= 0.05,]
results_13_p <- results_13[!is.na(results_13$padj) & results_13$padj<= 0.05,]
results_18_p <- results_18[!is.na(results_18$padj) & results_18$padj<= 0.05,]

summary(results_26)
sum(results_26$padj < 0.05, na.rm=TRUE)

## MLE used here, ask about MAP?
# Results log2 fold change (LFC) #
res_26_lfc <- lfcShrink(dds_deseq, coef = 2, res = results_26)
res_22_lfc <- lfcShrink(dds_deseq, coef = 2, res = results_22)


## MA-plot
plotMA(results_26, ylim = c(-2,2))
plotMA(res_26_lfc, ylim = c(-2,2))

plotMA(results_22, ylim = c(-2,2))
plotMA(res_22_lfc, ylim = c(-2,2))

rld <- rlog(dds, blind=FALSE) # or vts function 

# PCA plot
plotPCA(rld, intgroup=c("Compound"))

#pca_6 <- plotPCA(rld, intgroup=c("Compound"), returnData = TRUE)
# CFZ-22-t6-r3 seems to be an outlier

# Heat map
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## 24 hour ##
## ------- ##

dds_24 <-  DESeqDataSetFromMatrix(countData = CFZ_24hr,
                                  colData = coldata_24,
                                  design = ~ Compound)
dds_24

# Pre-filtering

# Keep only rows with reads over 10 counts, so that memory of object is reduced
keep <- rowSums(counts(dds_24)) >= 10
dds_24 <- dds_24[keep,]
dds_24

# Set the reference (control) level/ factor, against which others will be compared
dds_24$Compound <- relevel(dds_24$Compound, ref = "DMSO")

# Contrasting compounds against DMSO, using p adjusted value of 5% as cut off.
dds_deseq24 <- DESeq(dds_24)
results24_26 <- results(dds_deseq24, contrast = c("Compound", "26", "DMSO"), alpha = 0.05)
results24_26 <- results24_26[order(results24_26$padj),]
results24_22 <- results(dds_deseq24, contrast = c("Compound", "22", "DMSO"), alpha = 0.05)
results24_13 <- results(dds_deseq24, contrast = c("Compound", "13", "DMSO"), alpha = 0.05)
results24_18 <- results(dds_deseq24, contrast = c("Compound", "18", "DMSO"), alpha = 0.05)
results24_26_p <- results_26[!is.na(results24_26$padj) & results24_26$padj<= 0.05,]
results24_22_p <- results_22[!is.na(results24_22$padj) & results24_22$padj<= 0.05,]
results24_13_p <- results_13[!is.na(results24_13$padj) & results24_13$padj<= 0.05,]
results24_18_p <- results_18[!is.na(results24_18$padj) & results24_18$padj<= 0.05,]

summary(results24_26)
sum(results24_26$padj < 0.05, na.rm=TRUE)

## MLE used here, ask about MAP?
# Results log2 fold change (LFC) #
res24_26_lfc <- lfcShrink(dds_deseq24, coef = 2, res = results24_26)
res24_22_lfc <- lfcShrink(dds_deseq24, coef = 2, res = results24_22)


## MA-plot
plotMA(results24_26, ylim = c(-2,2))
plotMA(res24_26_lfc, ylim = c(-2,2))

plotMA(results24_22, ylim = c(-2,2))
plotMA(res24_22_lfc, ylim = c(-2,2))

rld24 <- rlog(dds_24, blind=FALSE) # or vts function 

# PCA plot
plotPCA(rld24, intgroup="Compound")
pca_24 <-plotPCA(rld24, intgroup=c("Compound"), returnData = TRUE)

# Heat map
sampleDists <- dist(t(assay(rld24)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



# test for getting PCA data
rv <- rowVars(assay(rld24))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]
pca <- prcomp(t(assay(rld24)[select, ]))

