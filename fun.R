library(DESeq2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gage)
library(gageData)

results_process <- function(dds, contrast, alpha, significant_only = TRUE, foldchange_threshold = FALSE){
  result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
  result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  # keep only non NA adjusted p values and <0.05
  if(significant_only == TRUE){
  result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  # All genes left now are 'statistically significant', so now order by biological significance 
  # Order by largest magnitude fold change
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
  }
  if(foldchange_threshold != FALSE){
    result_lfc <- result_lfc[abs(result_lfc$log2FoldChange) >= foldchange_threshold,]
    result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
  }
  final <- result_lfc
}

## Function for hashing gene id to name
gene_id_name <- function(ordered_results, top = 50, id_name_table = FALSE){
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  gene_keys <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
  #gene_keys <- gene_keys[!duplicated(gene_keys$ensembl_gene_id),]
  
  top_lfc <- ordered_results[1:top,]
  if(!all(rownames(top_lfc) %in% gene_keys$ensembl_gene_id)){
    print("not all genes in ensembl db")
    break
  }
  top_lfc$gene_name <- NA
  top_lfc$log_fold_change <- 0
  top_lfc$OverExpressed <- NA
  for(i in 1:nrow(top_lfc)){
    gene_id <- rownames(top_lfc)[i]
    name <- filter(gene_keys, gene_keys$ensembl_gene_id == gene_id)[1,2]
    top_lfc[i,ncol(top_lfc) -2] <- name
    # If positive log fold change, overexpressed = true
    top_lfc[i,ncol(top_lfc)] <- ifelse(top_lfc[i,2] >= 0,  TRUE, FALSE)
    top_lfc[i,ncol(top_lfc) -1] <- top_lfc[i,2]
  }
  # Output includes log fold change and p values etc. table
  final <- top_lfc[,-c(ncol(top_lfc),ncol(top_lfc)-1)]
  # output is just gene id vs gene name for the top statistically significant, largest log fold change genes
  if(id_name_table == TRUE){
    # top genes involved
    gene_top <- as.data.frame(cbind(rownames(top_lfc)[1:top], top_lfc$gene_name[1:top], top_lfc$log_fold_change[1:top] , top_lfc$OverExpressed[1:top]))
    colnames(gene_top) <- c("ensembl_gene_id", "hgnc_symbol", "log2FoldChange","OverExpressed")
    gene_top <- gene_top[!(is.na(gene_top$hgnc_symbol) | gene_top$hgnc_symbol==""), ]
    gene_top$log2FoldChange <- as.numeric(levels(gene_top$log2FoldChange))[gene_top$log2FoldChange]
   # rownames(gene_top) <- NULL
    final <- gene_top
  }
  return(final)
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


'%!in%' <- function(x,y)!('%in%'(x,y))

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

res_ouput <- function(result_lfc, gset, same_dir = TRUE, kegg_output = TRUE){
  result_lfc$gene_name <- convertIDs(row.names(result_lfc), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  result_lfc$entrez_id <- convertIDs(row.names(result_lfc), "ENSEMBL", "ENTREZID", org.Hs.eg.db)
  result_lfc <- result_lfc[!is.na(result_lfc$entrez_id) == TRUE,]
  if( kegg_output != TRUE){
    final <- as.data.frame(result_lfc[,c("gene_name", "entrez_id", "log2FoldChange")])
  }
  else{
  foldchanges <- result_lfc$log2FoldChange
  names(foldchanges) <- result_lfc$entrez_id
  keggres <- gage(foldchanges, gsets = gset, ref = NULL, samp = NULL, same.dir = same_dir)
  final <- list(keggres, foldchanges)
  }
  
  return(final)
}
dds_get <- function(gene_data, coldata, prefilter = 10, reference = "DMSO"){
  dds <-  DESeqDataSetFromMatrix(countData = gene_data,
                                 colData = coldata,
                                 design = ~ Compound)
  # Prefiltering, removing rows with sum of counts under 10 #
  keep <- rowSums(counts(dds)) >= prefilter
  dds <- dds[keep,]
  dds$Compound <- relevel(dds$Compound, ref = reference)
  
  dds_deseq <- DESeq(dds)
  return(dds_deseq)
}
pathway_full <- function(dds, contrast, gset, same_direction = TRUE, alpha = 0.05){
  # Getting latest human KEGG pathway
  result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
  result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)

  result_lfc$gene_name <- convertIDs(row.names(result_lfc), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  result_lfc$entrez_id <- convertIDs(row.names(result_lfc), "ENSEMBL", "ENTREZID", org.Hs.eg.db)
  result_lfc <- result_lfc[!is.na(result_lfc$entrez_id) == TRUE,]
  foldchanges <- result_lfc$log2FoldChange
  names(foldchanges) <- result_lfc$entrez_id
  keggres <- gage(foldchanges, gsets = gset, ref = NULL, samp = NULL, same.dir = same_direction)
  final <- list(keggres, foldchanges)
  
  return(final)
  
}