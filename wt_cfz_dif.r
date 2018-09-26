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

# 6 hours first, no outliers
sig_cfz_26_6_up <- res_26_6_useful[!is.na(res_26_6_useful$padj) & res_26_6_useful$padj < 0.05 & res_26_6_useful$log2FoldChange >1,]
sig_wt_26_6_up <- res_26_6_useful_wt[!is.na(res_26_6_useful_wt$padj) & res_26_6_useful_wt$padj < 0.05 & res_26_6_useful_wt$log2FoldChange >1,]

sig_cfz_26_6_down <- res_26_6_useful[!is.na(res_26_6_useful$padj) & res_26_6_useful$padj < 0.05 & res_26_6_useful$log2FoldChange < -1,]
sig_wt_26_6_down <- res_26_6_useful_wt[!is.na(res_26_6_useful_wt$padj) & res_26_6_useful_wt$padj < 0.05 & res_26_6_useful_wt$log2FoldChange < -1,]


sig_cfz_13_6_up <- res_13_6_useful[!is.na(res_13_6_useful$padj) & res_13_6_useful$padj < 0.05 & res_13_6_useful$log2FoldChange >1,]
sig_wt_13_6_up <- res_13_6_useful_wt[!is.na(res_13_6_useful_wt$padj) & res_13_6_useful_wt$padj < 0.05 & res_13_6_useful_wt$log2FoldChange >1,]

sig_cfz_13_6_down <- res_13_6_useful[!is.na(res_13_6_useful$padj) & res_13_6_useful$padj < 0.05 & res_13_6_useful$log2FoldChange < -1,]
sig_wt_13_6_down <- res_13_6_useful_wt[!is.na(res_13_6_useful_wt$padj) & res_13_6_useful_wt$padj < 0.05 & res_13_6_useful_wt$log2FoldChange < -1,]


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


## Top 50 genes for each

res_26_6_useful_wt$OverExpressed <- ifelse(res_26_6_useful_wt$log2FoldChange >=0, TRUE, FALSE)
res_26_6_useful_wt <- as.data.frame(res_26_6_useful_wt)
res_26_6_useful_wt$OverExpressed <- as.factor(res_26_6_useful_wt$OverExpressed)
g<-ggplot(data=res_26_6_useful_wt[1:50,], aes(x=gene_name, y=log2FoldChange, fill = OverExpressed)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(size=rel(1.3), angle=90)) + xlab("Gene") + ylab("Log2 Fold Change") +
  scale_fill_manual(values = c("Red", "Black")) 
print(g)


samples_list <- list(res_26_6_useful_wt, res_26_24_useful_wt, res_13_6_useful_wt, res_13_24_useful_wt)
gene_list <- unique(c(res_26_6_useful_wt$gene_name[1:15],res_26_24_useful_wt$gene_name[1:15],res_13_6_useful_wt$gene_name[1:15], res_13_24_useful_wt$gene_name[1:15]))
plot_table <- data.frame()

plot_table[1:2,1] <- "26"
plot_table[3:4,1] <- "13"
plot_table[c(1,3),2] <- "6"
plot_table[c(2,4),2] <- "24"
rownames(plot_table) <- c("26_6", "26_24","13_6", "13_24" )
colnames(plot_table) <- c("Compound", "TimePoint")

for(i in seq_along(samples_list)){
  for(j in seq_along(gene_list)){
    subset_data <- samples_list[[i]] %>% filter(str_detect(gene_name, gene_list[j]))
    plot_table[i,j + 2] <- subset_data[1,"log2FoldChange"]
  }
}
colnames(plot_table) <- c("Compound", "TimePoint",gene_list)
plot_table$Compound <- as.factor(plot_table$Compound)
plot_table$TimePoint <- as.factor(plot_table$TimePoint)

plot_table_tidy <- gather(plot_table,"gene", "Log2FoldChange", -c(Compound, TimePoint))

g <- ggplot(data=plot_table_tidy, aes(x=gene, y=Log2FoldChange, fill = Compound)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge())+ theme(axis.text.x = element_text(size=rel(1), angle=90))
print(g)


## Joint plots

CFZ_genes <- read.csv("CFZ_genes.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("coldata.csv", row.names = 1)

# outlier 24hr: CFZ-DMSO-t24-r2, CFZ-18-t24-r3
outlier_row_column_list <- c(which(rownames(coldata) == "CFZ-18-t24-r3"),which(rownames(coldata) == "CFZ-DMSO-t24-r2"))
CFZ_genes <- CFZ_genes[, -outlier_row_column_list]
coldata_compare <- coldata[-outlier_row_column_list,]
coldata_compare$CellType <- "CFZ"

# Load WT
WT_genes <- read.csv("WT_genes.csv", row.names = 1, check.names = FALSE)
WT_genes_no_carf <- dplyr::select(WT_genes, -contains("carf"))
coldata_wt <- read.csv("coldata_wt.csv", row.names = 1)
coldata_wt$Compound <- as.factor(coldata_wt$Compound)
coldata_wt$CellType <- "WT"
coldata_wt_compare <- coldata_wt[coldata_wt$Compound != "carf",]

# Bind
all_genes <- cbind(WT_genes_no_carf,CFZ_genes)
coldata_all <- rbind(coldata_wt_compare, coldata_compare)
coldata_all$CellType <- as.factor(coldata_all$CellType)
coldata_all$TimePoint <- as.factor(coldata_all$TimePoint)

all(rownames(coldata_all) == colnames(all_genes)) # TRUE

# DESeq2
dds_compare_all <-  DESeqDataSetFromMatrix(countData = all_genes,
                                   colData = coldata_all,
                                   design = ~ TimePoint + CellType + Compound)
# Prefiltering, removing rows with sum of counts under 10 #
keep <- rowSums(counts(dds_compare_all)) >= 10
dds_compare_all <- dds_compare_all[keep,]
dds_compare_all$CellType <- relevel(dds_compare_all$CellType, ref = "WT")
dds_compare_all$Compound <- relevel(dds_compare_all$Compound, ref = "DMSO")
dds_compare_all$TimePoint <- relevel(dds_compare_all$TimePoint, ref = "6")

dds_compare_all <- DESeq(dds_compare_all)
vst_all <- vst(dds_compare_all, blind=FALSE)

plotPCA(vst_all, intgroup=c("CellType","Compound"))

