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
