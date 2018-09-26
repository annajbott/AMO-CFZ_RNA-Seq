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

## WT vs CFZ DMSO ##
## -------------- ##

# Motivation to check drug resistance #


# Load WT
WT_genes <- read.csv("WT_genes.csv", row.names = 1, check.names = FALSE)
coldata_wt <- read.csv("coldata_wt.csv", row.names = 1)
coldata_wt$Compound <- as.factor(coldata_wt$Compound)

## Separate time points ##
## -------------------- ##

WT_6hr <- dplyr::select(WT_genes,-contains("t24"))
WT_24hr <- dplyr::select(WT_genes,-contains("t6"))
coldata_wt_6hr <- coldata_wt[coldata_wt$TimePoint == 6,]
coldata_wt_24hr <- coldata_wt[coldata_wt$TimePoint == 24,]

# Check in same order
all(colnames(WT_6hr) == rownames(coldata_wt_6hr)) # True, same order
all(colnames(WT_24hr) == rownames(coldata_wt_24hr)) # True, same order

## DESeq2 ##
## ------ ##

dds_deseq_wt <- dds_get(gene_data = WT_6hr, coldata = coldata_wt_6hr, prefilter = 10, reference = "DMSO")
dds_deseq24_wt <- dds_get(gene_data = WT_24hr, coldata = coldata_wt_24hr, prefilter = 10, reference = "DMSO")

res_26_lfc_wt <- results_process(dds_deseq_wt, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_22_lfc_wt <- results_process(dds_deseq_wt, contrast = c("Compound", "22", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_13_lfc_wt <- results_process(dds_deseq_wt, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_18_lfc_wt <- results_process(dds_deseq_wt, contrast = c("Compound", "18", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_carf_lfc_wt <- results_process(dds_deseq_wt, contrast = c("Compound", "carf", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )


# Taken out a replicate so padj values are using sd with only two replicates, instead set threshold as 1 log2 fold change (double expression)
res_26_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_22_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "22", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_13_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_18_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "18", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_carf_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "carf", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )


res_26_6_useful_wt <- res_useful(res_26_lfc_wt, tf = TRUE)
res_22_6_useful_wt <- res_useful(res_22_lfc_wt, tf = TRUE)
res_13_6_useful_wt <- res_useful(res_13_lfc_wt, tf = TRUE)
res_18_6_useful_wt <- res_useful(res_18_lfc_wt, tf = TRUE)
res_carf_6_useful_wt <- res_useful(res_carf_lfc_wt, tf = NULL)

res_26_24_useful_wt <- res_useful(res_26_lfc_24_wt, tf = TRUE)
res_22_24_useful_wt <- res_useful(res_22_lfc_24_wt, tf = TRUE)
res_13_24_useful_wt <- res_useful(res_13_lfc_24_wt, tf = TRUE)
res_18_24_useful_wt <- res_useful(res_18_lfc_24_wt, tf = TRUE)
res_carf_24_useful_wt <- res_useful(res_carf_lfc_24_wt, tf = NULL)


# Data object so don't have to run the above
save(dds_deseq_wt, dds_deseq24_wt, res_26_6_useful_wt, res_22_6_useful_wt, res_13_6_useful_wt, res_18_6_useful_wt,res_26_24_useful_wt ,res_22_24_useful_wt , res_13_24_useful_wt , res_18_24_useful_wt, file = "WT_useful_deseq_results.RData")

# Can run from here instead- shortcut
load("WT_useful_deseq_results.RData")

eTerm_dmso_wt_26_6 <- enricher_analysis(dds_deseq_wt, c("Compound", "26", "DMSO"), ontology = "MsigdbC2REACTOME", result_lfc = res_26_6_useful_wt, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500)
eTerm_dmso_wt_13_6 <- enricher_analysis(dds_deseq_wt, c("Compound", "13", "DMSO"), ontology = "MsigdbC2REACTOME", result_lfc = res_13_6_useful_wt, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500)
eTerm_dmso_wt_26_24 <- enricher_analysis(dds_deseq24_wt, c("Compound", "26", "DMSO"), ontology = "MsigdbC2REACTOME", result_lfc = res_26_24_useful_wt, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500)
eTerm_dmso_wt_13_24 <- enricher_analysis(dds_deseq24_wt, c("Compound", "13", "DMSO"), ontology = "MsigdbC2REACTOME", result_lfc = res_13_24_useful_wt, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500)


subnet_26_6_wt <- subneter_analysis(dds_deseq, c("Compound", "26", "DMSO"), result_lfc = res_26_6_useful_wt, alpha = 0.05, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_26_6_wt, pattern=-log10(as.numeric(V(subnet_26_6_wt)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.5, vertex.label.cex = 0.6, vertex.label.color = "black")
#
subnet_26_24_wt <- subneter_analysis(dds_deseq24, c("Compound", "26", "DMSO"), result_lfc = res_26_24_useful_wt, alpha = 0.05, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_26_24_wt, pattern=-log10(as.numeric(V(subnet_26_24_wt)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
subnet_13_6_wt <- subneter_analysis(dds_deseq, c("Compound", "13", "DMSO"), result_lfc = res_13_6_useful_wt, alpha = 0.05, number_top_genes = 500,  subnet.size = 50)
xVisNet(g=subnet_13_6_wt, pattern=-log10(as.numeric(V(subnet_13_6_wt)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
subnet_13_24_wt <- subneter_analysis(dds_deseq24, c("Compound", "13", "DMSO"), result_lfc = res_13_24_useful_wt, alpha = 0.05, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_13_24_wt, pattern=-log10(as.numeric(V(subnet_13_24_wt)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#

