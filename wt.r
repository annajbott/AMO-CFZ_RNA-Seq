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

# Taken out a replicate so padj values are using sd with only two replicates, instead set threshold as 1 log2 fold change (double expression)
res_26_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_22_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "22", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE)
res_13_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )
res_18_lfc_24_wt <- results_process(dds_deseq24_wt, contrast = c("Compound", "18", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = FALSE )

res_26_6_useful_wt <- res_useful(res_26_lfc_wt, tf = TRUE)
res_22_6_useful_wt <- res_useful(res_22_lfc_wt, tf = TRUE)
res_13_6_useful_wt <- res_useful(res_13_lfc_wt, tf = TRUE)
res_18_6_useful_wt <- res_useful(res_18_lfc_wt, tf = TRUE)
res_26_24_useful_wt <- res_useful(res_26_lfc_24_wt, tf = TRUE)
res_22_24_useful_wt <- res_useful(res_22_lfc_24_wt, tf = TRUE)
res_13_24_useful_wt <- res_useful(res_13_lfc_24_wt, tf = TRUE)
res_18_24_useful_wt <- res_useful(res_18_lfc_24_wt, tf = TRUE)

# Data object so don't have to run the above
save(res_26_6_useful_wt, res_22_6_useful_wt, res_13_6_useful_wt, res_18_6_useful_wt,res_26_24_useful_wt ,res_22_24_useful_wt , res_13_24_useful_wt , res_18_24_useful_wt, file = "WT_useful_deseq_results.RData")

# Can run from here instead
load("WT_useful_deseq_results.RData")