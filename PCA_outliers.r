library(tidyverse)
library("DESeq2")
library("RColorBrewer")
library(pheatmap)
library("IHW")
library(biomaRt)
source("fun.R")
library(gage)
library(gageData)
library(pathview)
library(org.Hs.eg.db)


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
outlier_row_column_list <- c(which(rownames(coldata) == "CFZ-18-t24-r3"),which(rownames(coldata) == "CFZ-DMSO-t24-r2"))
CFZ_genes <- CFZ_genes[, -outlier_row_column_list]
coldata <- coldata[-outlier_row_column_list,]

# individual time point, subset of above data

CFZ_6hr <- CFZ_genes[,1:15]
# NCP 22 r3 isn't really an outlier, onlyon PC2, which accounts for 8% variance
coldata_6 <- coldata[1:15,]
CFZ_24hr <- CFZ_genes[,16:28]
coldata_24 <- coldata[16:28,]

## Using DESeq2 Package ##
## -------------------- ##

dds_deseq <- dds_get(gene_data = CFZ_6hr, coldata = coldata_6, prefilter = 10, reference = "DMSO")
dds_deseq24 <- dds_get(gene_data = CFZ_24hr, coldata = coldata_24, prefilter = 10, reference = "DMSO")

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

## Differential expression analysis without outliers ##
## ------------------------------------------------- ##

res_26_lfc <- results_process(dds_deseq, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_22_lfc <- results_process(dds_deseq, contrast = c("Compound", "22", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_13_lfc <- results_process(dds_deseq, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_18_lfc <- results_process(dds_deseq, contrast = c("Compound", "18", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )

# Taken out a replicate so padj values are using sd with only two replicates, instead set threshold as 1 log2 fold change (double expression)
res_26_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_22_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "22", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_13_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_18_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "18", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )

######################
## Pathway Analysis ##
######################

# Might be needed later???
out.suffix <- "deseq2"

## Using GATE and KEGG ##
## ------------------- ##

# Getting latest human KEGG pathway
kg.hsa <- kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)
kegg.sigmet <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]

# List of results, followed by fold change data frame
kegg_26_6 <- pathway_full(dds = dds_deseq, contrast = c("Compound", "26", "DMSO"), gset = kegg.sigmet,  same_direction = TRUE)
kegg_26_24 <- pathway_full(dds = dds_deseq24, contrast = c("Compound", "26", "DMSO"), gset = kegg.sigmet,  same_direction = TRUE)
kegg_13_6 <-  pathway_full(dds = dds_deseq, contrast = c("Compound", "13", "DMSO"), gset = kegg.sigmet,  same_direction = TRUE)
kegg_13_24 <- pathway_full(dds = dds_deseq24, contrast = c("Compound", "13", "DMSO"), gset = kegg.sigmet,  same_direction = TRUE)
# Test
lapply(kegg_13_24[[1]], head)

kegg_26 <- kegg_26_6[[1]]
fc_kegg_26 <- kegg_26_6[[2]]


gs <- unique(unlist(kegg.sigmet[rownames(kegg_26$greater)[1]]))
fc_kegg_26 <- cbind(fc_kegg_26) # From biostar forum
essData <- essGene(gs, fc_kegg_26, ref =NULL, samp =NULL)
convertIDs(rownames(as.data.frame(essData)),"ENTREZID","SYMBOL", org.Hs.eg.db)

# this is for heat map probably won't work with deseq2 format
head(essData, 4)
ref1=1:6
samp1=7:12
#generated text file for data table, pdf files for heatmap and scatterplot
for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
           samp = samp1, outname = outname, txt = TRUE, heatmap = TRUE,
           Colv = FALSE, Rowv = FALSE, dendrogram = "none", limit = 3, scatterplot = TRUE)
}



## Pathview crap
# Get the pathways
kegg_pathways_26_6 = data.frame(id=rownames(kegg_26_6[[1]]$greater), kegg_26_6[[1]]$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=3) %>% 
  .$id %>% 
  as.character()
# Get the IDs.
kegg_ids_26_6 = substr(kegg_pathways_26_6, start=1, stop=8)
# 13_6 dpwn_regulated
kegg_pathways_13_24 = data.frame(id=rownames(kegg_13_24[[1]]$less), kegg_13_24[[1]]$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=1) %>% 
  .$id %>% 
  as.character()
# Get the IDs.
kegg_ids_13_24 = substr(kegg_pathways_13_24, start=1, stop=8)

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(kegg_ids_13_24, function(pid) pathview(gene.data=kegg_13_24[[2]], pathway.id=pid, species="hsa", out.suffix=out.suffix))

