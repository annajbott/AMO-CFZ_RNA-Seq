library(tidyverse)
source("fun.R")
library(org.Hs.eg.db)
library(stringr)

####################################################
## Compare WT dif expressed and CFZ dif expressed ##
####################################################

# Load stuff from data objects (no long waiting)

load("WT_useful_deseq_results.RData")
load("useful_deseq_results.RData")
# Everything, can filter by padj and fold change or whether TF 
# Format: res_compound_time_useful or res_compound_time_useful_wt

# Check specified genes (example)
res_26_6_useful %>% filter(str_detect(gene_name, "HIST"))

# Check how many differentially expressed genes and percent upregulated (6 and 24hrs)
up_regulated(res_18_6_useful_wt,res_18_24_useful_wt)

res_26_6_useful_wt[!is.na(res_26_6_useful_wt$padj) & res_26_6_useful_wt$padj < 0.05 & abs(res_26_6_useful_wt$log2FoldChange) >1,]
res_26_6_useful[!is.na(res_26_6_useful$padj) & res_26_6_useful$padj < 0.05 & abs(res_26_6_useful$log2FoldChange) >1,]
#
res_26_24_useful_wt[!is.na(res_26_6_useful_wt$padj) & res_26_6_useful_wt$padj < 0.05 & abs(res_26_6_useful_wt$log2FoldChange) >1,]
setdiff(res_26_24_useful[abs(res_26_24_useful$log2FoldChange) >1 & res_26_24_useful$tf == TRUE,]$gene_name,res_26_24_useful_wt[abs(res_26_24_useful_wt$log2FoldChange) >1 & res_26_24_useful_wt$tf == TRUE,]$gene_name)

setdiff(res_26_24_useful_wt[abs(res_26_24_useful_wt$log2FoldChange) >1 & res_26_24_useful_wt$tf == TRUE,]$gene_name,res_26_24_useful[abs(res_26_24_useful$log2FoldChange) >0.5 & res_26_24_useful$tf == TRUE,]$gene_name)
