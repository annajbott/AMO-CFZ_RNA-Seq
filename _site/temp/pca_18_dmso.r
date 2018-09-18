library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")
source("fun.R")
library(stringr)
library(ggplot2)

# 18 and 24 variation

CFZ_genes <- read.csv("CFZ_genes.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("coldata.csv", row.names = 1)

# From PCA plots and data from 'returnData' pca option
# outlier 6hr: CFZ-22-t6-r3
# outlier 24hr: CFZ-DMSO-t24-r2, CFZ-18-t24-r3, ?CFZ-26-t24-r1?
outlier_row_column_list <- c(which(rownames(coldata) == "CFZ-26-t24-r1"))
CFZ_genes <- CFZ_genes[, -outlier_row_column_list]
coldata <- coldata[-outlier_row_column_list,]

CFZ_24hr <- CFZ_genes[,16:29]
coldata_24 <- coldata[16:29,]

dds_24_18_dmso <-  DESeqDataSetFromMatrix(countData = CFZ_24hr,
                                  colData = coldata_24,
                                  design = ~ Compound)
dds_24_18_dmso

## Pre-filtering ##
# Keep only rows with reads over 10 counts, so that memory of object is reduced
keep <- rowSums(counts(dds_24_18_dmso)) >= 10
dds_24_18_dmso <- dds_24_18_dmso[keep,]
dds_24_18_dmso

# Set the reference (control) level/ factor, against which others will be compared
dds_24_18_dmso$Compound <- relevel(dds_24_18_dmso$Compound, ref = "DMSO")

# Contrasting compounds against DMSO, using p adjusted value of 5% as cut off.
dds_24_18_dmso <- DESeq(dds_24_18_dmso)


rld24_18_dmso <- rlog(dds_24_18_dmso, blind=FALSE) 

plotPCA(rld24_18_dmso, intgroup="Compound")

# 2 outliers from 24 hour PCA, 1 DMSO and 1 MAZ 18. Affects PC1
rv <- rowVars(assay(rld24_18_dmso))
select_rv <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                      length(rv)))]
pca <- prcomp(t(assay(rld24_18_dmso)[select_rv, ]))

pca_rot <- as.data.frame(pca$rotation[,1:2])
PC1_18_dmso <- pca_rot[rev(order(pca_rot$PC1)),][1:100,] # Order by genes contributing most variation to PC1
symbols <- convertIDs(rownames(PC1_18_dmso), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
PC1_18_dmso$symbols <- convertIDs(rownames(PC1_18_dmso), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
PC1_18_dmso <- PC1_18_dmso[,c("symbols","PC1")]
PC1_18_dmso <- PC1_18_dmso[!is.na(PC1_18_dmso$symbols),]

eTerm_PC1 <- xEnricherGenes(data = PC1_18_dmso$symbols, ontology = "MsigdbC2KEGG")
eTerm_PC1_re <- xEnricherGenes(data = PC1_18_dmso$symbols, ontology = "REACTOME")
#lapply(eTerm_PC1_re$annotation,head)$
  
list_eTerm <- list(eTerm_PC1, eTerm_PC1_re)
names(list_eTerm) <- c('KEGG','REACTOME')
bp_Pathway <- xEnrichCompare(list_eTerm, displayBy="fc", FDR.cutoff=5e-2, wrap.width=50)
bp_Pathway + theme(axis.text.y=element_text(size=10))
  