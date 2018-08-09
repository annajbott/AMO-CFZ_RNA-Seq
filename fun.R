library(DESeq2)

results_process <- function(dds, contrast, alpha){
  result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
  result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  # keep only non NA adjusted p values and <0.05
  result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  # Order by largest magnitude fold change
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
}