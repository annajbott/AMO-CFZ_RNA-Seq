library(tidyverse)
library("DESeq2")
library(readxl)
library("RColorBrewer")
library(pheatmap)


gene_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/genes.tsv.gz")
gene_tib2 <- as.data.frame(gene_tib)[-1]
rownames(gene_tib2) <- as.data.frame(gene_tib)[,1]
head(gene_tib2)

#transcripts_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/transcripts.tsv.gz")

# Looking at data structure
glimpse(gene_tib2)
summary(gene_tib2)

# Remove sample number at the end
colnames(gene_tib2) <- sub("_[^_]+$", "", colnames(gene_tib2))
# Remove index number
colnames(gene_tib2) <-gsub("^i[0-9][0-9]-","",colnames(gene_tib2))


# Genes going down, sample name going along
CFZ_genes <- select(gene_tib2, -contains("WT"))


# Separate by time points, might need to move later
CFZ_6hr <- select(CFZ_genes, -contains("t24"))
CFZ_24hr <- select(CFZ_genes, -contains("t6"))

# Read in sample sheet
#coldata <- read_excel("/home/osboxes/Documents/analysis/AMO-CFZ_RNA-Seq/sample_sheet.xlsx")
#names(coldata)[1] <- "SampleName"
#coldata <- read.csv("/home/osboxes/Documents/analysis/AMO-CFZ_RNA-Seq/sample_sheet.xlsx", row.names =1)

coldata <- read.csv("/home/osboxes/Documents/analysis/AMO-CFZ_RNA-Seq/sample_sheet2.csv", row.names =1)
coldata <- coldata[-1]
names(coldata) <- c("Compound", "TimePoint")

# Replace underscores for hyphens
rownames(coldata) <- gsub("_","-",rownames(coldata))


# Ensure the rows of sample sheet match the columns of the gene matrix
rownames(coldata)
colnames(CFZ_genes)
all(rownames(coldata) %in% colnames(CFZ_genes)) # True 

# They have the same values but not the same order
all(rownames(coldata) == colnames(CFZ_genes)) # False
# Put them in the same order as the sample sheet
CFZ_genes <- CFZ_genes[, rownames(coldata)]
# Test they are in the same order now
all(rownames(coldata) == colnames(CFZ_genes)) # Now true

## Separate into 6 and 24 hour time points ##
## --------------------------------------- ##

CFZ_6hr <- select(CFZ_genes, -contains("t24"))
CFZ_24hr <- select(CFZ_genes, -contains("t6"))

## 6 and 24 hour versions of smaple sheet, keep row names
coldata_6<-coldata[coldata$TimePoint == 6, ]
coldata_24 <- coldata[coldata$TimePoint == 24, ]

# Check 6 and 24 hour filtered tables also match
all(rownames(coldata_6) == colnames(CFZ_6hr)) # TRUE
all(rownames(coldata_24) == colnames(CFZ_24hr)) # TRUE

## DESeq2 data set ##
## --------------- ##
# 6hr
dds <-  DESeqDataSetFromMatrix(countData = CFZ_6hr,
                               colData = coldata_6,
                               design = ~ Compound)
dds

###################
## Pre-filtering ##
###################
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

rld <- rlog(dds, blind=FALSE)

plotPCA(rld, intgroup=c("Compound"))

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
