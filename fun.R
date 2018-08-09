library(DESeq2)
library(biomaRt)

results_process <- function(dds, contrast, alpha){
  result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
  result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  # keep only non NA adjusted p values and <0.05
  result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  # All genes left now are 'statistically significant', so now order by biological significance 
  # Order by largest magnitude fold change
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
}

## Function for hashing gene id to name
gene_id_name <- function(ordered_results, top = 50){
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  gene_keys <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
  #gene_keys <- gene_keys[!duplicated(gene_keys$ensembl_gene_id),]
  
  top_lfc <- ordered_results[1:top,]
  if(!all(rownames(top_lfc) %in% gene_keys$ensembl_gene_id)){
    print("not all genes in ensembl db")
    break
  }
  top_lfc$gene_name <- NA
  for(i in 1:nrow(top_lfc)){
    gene_id <- rownames(top_lfc)[i]
    name <- filter(gene_keys, gene_keys$ensembl_gene_id == gene_id)[1,2]
    top_lfc[i,ncol(top_lfc)] <- name
  }
  return(top_lfc)
}

gene_id_name_raw <- function(string_vector){
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  gene_keys <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
  data_frame_genes <- data.frame(matrix(NA, nrow = length(string_vector), ncol = 2))
  for( i in 1:length(string_vector)){
    data_frame_genes[i,1] <- string_vector[i]
    data_frame_genes[i,2] <- filter(gene_keys, gene_keys$ensembl_gene_id == string_vector[i])[1,2]
  }
  colnames(data_frame_genes) <- c("ensembl_gene_id", "hgnc_symbol")
  data_frame_genes <- data_frame_genes[!(is.na(data_frame_genes$hgnc_symbol) | data_frame_genes$hgnc_symbol==""), ]
  return(data_frame_genes)
}