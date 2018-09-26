library(DESeq2)
#library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
#library(gage)
#library(gageData)
library("XGR")
library(tidyverse)

results_process <- function(dds, contrast, alpha, significant_only = TRUE, foldchange_threshold = FALSE, order = TRUE){
  result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
  result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  # keep only non NA adjusted p values and <0.05
  if(significant_only == TRUE){
  result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  # All genes left now are 'statistically significant', so now order by biological significance 
  }
  if(foldchange_threshold != FALSE){
    result_lfc <- result_lfc[abs(result_lfc$log2FoldChange) >= foldchange_threshold,]
  }
  if(order == TRUE){
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

res_ouput <- function(result_lfc, gset, same_dir = TRUE, kegg_output = TRUE, alpha = 1){
  result_lfc$gene_name <- convertIDs(row.names(result_lfc), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  result_lfc$entrez_id <- convertIDs(row.names(result_lfc), "ENSEMBL", "ENTREZID", org.Hs.eg.db)
  result_lfc <- result_lfc[!is.na(result_lfc$entrez_id) == TRUE & !is.na(result_lfc$padj)& result_lfc$padj<= alpha,]
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
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
dds_get <- function(gene_data, coldata, prefilter = 10, reference = "DMSO", cell_type = FALSE){
  if(cell_type == TRUE ){
    dds <-  DESeqDataSetFromMatrix(countData = gene_data,
                                   colData = coldata,
                                   design = ~CellType)
    keep <- rowSums(counts(dds)) >= prefilter
    dds <- dds[keep,]
    dds$CellType <- relevel(dds$CellType, ref = reference)
  }else{
    dds <-  DESeqDataSetFromMatrix(countData = gene_data,
                                   colData = coldata,
                                   design = ~Compound)
    # Prefiltering, removing rows with sum of counts under 10 #
    keep <- rowSums(counts(dds)) >= prefilter
    dds <- dds[keep,]
    dds$Compound <- relevel(dds$Compound, ref = reference)
  }

  
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
# For XGR pathway analysis
enricher_analysis <- function(dds, contrast, ontology, result_lfc = NULL, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200, same.dir = TRUE){
  if(is.null(result_lfc)){
    result_noshrink <- results(dds, contrast = contrast, alpha = alpha)
    result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  }
  background_symbols <- convertIDs(rownames(result_lfc),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
  background_symbols <-  as.vector(background_symbols[!is.na(background_symbols)])
  
  if(foldchange_threshold == TRUE){
    result_lfc <- result_lfc[abs(result_lfc$log2FoldChange) >= foldchange_threshold,]
  }
  else{
    result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  }
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
  if(same.dir != TRUE){
    # Separate up and down regulated genes
    result_lfc_up <- result_lfc[result_lfc$log2FoldChange > 0,]
    result_lfc_down <- result_lfc[result_lfc$log2FoldChange < 0,]
    ifelse(number_top_genes > length(rownames(result_lfc_up)),number_top_genes_up <- length(rownames(result_lfc_up)), number_top_genes_up <- number_top_genes)
    ifelse(number_top_genes > length(rownames(result_lfc_down)), number_top_genes_down <- length(rownames(result_lfc_down)), number_top_genes_down <- number_top_genes)
    data_symbols_up <- convertIDs(rownames(result_lfc_up)[1:number_top_genes_up],"ENSEMBL", "SYMBOL", org.Hs.eg.db)
    data_symbols_up <- as.vector(data_symbols_up[!is.na(data_symbols_up)])
    data_symbols_down <- convertIDs(rownames(result_lfc_down)[1:number_top_genes_down],"ENSEMBL", "SYMBOL", org.Hs.eg.db)
    data_symbols_down <- as.vector(data_symbols_down[!is.na(data_symbols_down)])
    eTerm_up <- xEnricherGenes(data = data_symbols_up, background = background_symbols, ontology = ontology)
    eTerm_down <- xEnricherGenes(data = data_symbols_down, background = background_symbols, ontology = ontology)
    final <- list(xEnrichConciser(eTerm_up), xEnrichConciser(eTerm_down))
  }else{
  # All genes, up or down regulated
  if(number_top_genes > length(rownames(result_lfc))) {number_top_genes <- length(rownames(result_lfc))}
  data_symbols <- convertIDs(rownames(result_lfc)[1:number_top_genes],"ENSEMBL", "SYMBOL", org.Hs.eg.db)
  data_symbols <- as.vector(data_symbols[!is.na(data_symbols)])
  eTerm <- xEnricherGenes(data = data_symbols, background = background_symbols, ontology = ontology)
  final <- xEnrichConciser(eTerm, cutoff = c(0.9, 0.5), verbose = F)
  }
  return(final)
}


# For XGR network subnetergenes
subneter_analysis <- function(dds, contrast, result_lfc = NULL, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500, network = "STRING_high", subnet.size=NULL){
  if(is.null(result_lfc)){
    result_noshrink <- results(dds, contrast = contrast, alpha = 0.05)
    result_lfc <- lfcShrink(dds, contrast = contrast, res = result_noshrink)
  }
  # Use FC instead of p values
  if(foldchange_threshold == TRUE){
    result_lfc <- result_lfc[abs(result_lfc$log2FoldChange) >= foldchange_threshold,]
  }
  else if(alpha != 1){
    result_lfc <- result_lfc[!is.na(result_lfc$padj) & result_lfc$padj<= alpha,]
  }
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]

  
  
  result_lfc$data_symbols <- convertIDs(rownames(result_lfc),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
  result_lfc <- result_lfc[!is.na(result_lfc$data_symbols),]
  if(number_top_genes > length(rownames(result_lfc))) {number_top_genes <- length(rownames(result_lfc))}
  data <- as.data.frame(result_lfc[,c("data_symbols", "padj")])
  rownames(data) <- NULL
  data <- data[1:number_top_genes,]
  # Significance threshold given as default 0.01, if subnet size given will overwrite default of NULL to size and assign new threshold
  subnet <- xSubneterGenes(data= data, network= network, subnet.significance = 0.01, subnet.size= subnet.size)
  return(subnet)
}

find_tf <- function(gene_list, tf_reference){
  tf_gene_list <- c()
  for(i in gene_list){
    if(any(tf_reference == i )){
      tf_gene_list <- c(tf_gene_list,i)
    }
  }
  return(tf_gene_list)
}

res_useful <- function(result_lfc, alpha = 1, tf = NULL){
  result_lfc <- result_lfc[,c(1,2,6)]
  result_lfc$gene_name <- convertIDs(row.names(result_lfc), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  result_lfc$entrez_id <- convertIDs(row.names(result_lfc), "ENSEMBL", "ENTREZID", org.Hs.eg.db)
  result_lfc <- result_lfc[!is.na(result_lfc$entrez_id) == TRUE ,]
  if(alpha < 1){
    result_lfc <- result_lfc[!is.na(result_lfc$padj)& result_lfc$padj<= alpha,]
  }
  result_lfc <- result_lfc[rev(order(abs(result_lfc$log2FoldChange))),]
  if(!is.null(tf)){
    TF_csv <- read_csv("temp/TFCheckpoint.csv")
    hs_tf <- dplyr::select(TF_csv, c("gene_symbol", "entrez_human", "gene_name", "synonym", "DbTF"))
    result_lfc$tf <- NULL
    for(i in 1:length(rownames(result_lfc))){
      result_lfc$tf[i] <- any(hs_tf$gene_symbol == result_lfc$gene_name[i] )
    }
    result_lfc <- as.data.frame(result_lfc)
    result_lfc <- result_lfc[,c(4,5,1,2,3,6)]
  }else{    
    result_lfc <- as.data.frame(result_lfc)
    result_lfc <- result_lfc[,c(4,5,1,2,3)]}

  return(result_lfc)
}

up_regulated <- function(df1, df2){
  
  df1_sig <- df1[!is.na(df1$padj) & df1$padj < 0.05 ,]
  df2_sig <- df2[!is.na(df2$padj) & df2$padj < 0.05 ,]
  x <- c(df1_sig$gene_name, df2_sig$gene_name)
  total <- n_distinct(x)
  up1 <- n_distinct(df1_sig[df1_sig$log2FoldChange >0 ,])
  down1 <- n_distinct(df1_sig[df1_sig$log2FoldChange <0 ,])
  percent1 <- up1/(up1+down1)
  up2 <- n_distinct(df2_sig[df2_sig$log2FoldChange >0 ,])
  down2 <- n_distinct(df2_sig[df2_sig$log2FoldChange <0 ,])
  percent2 <- up2/(up2+down2)
  average_percent <- (percent1 + percent2)/2
  if(length(rownames(df2_sig)) == 0 & length(rownames(df1_sig)) == 0){
    average_percent <- NULL
  }else if(length(rownames(df2_sig)) == 0){
    average_percent <- percent1
  }else if(length(rownames(df1_sig)) == 0){
    average_percent <- percent2
  }
  return(list(total,average_percent))
}

enriched_genes_annotation <- function(eterm, pathway_entrez_array, results_useful, sig_only = FALSE){
  ## TCR signalling
  pathway_name_genes <- c()
  pathway_name_fc <- c()
  for(i in 1:length(pathway_entrez_array)){
    row <- res_dmso_624_useful[res_dmso_624_useful$entrez_id == pathway_entrez_array[[i]],]
    if(length(rownames(row)) == 1){
      pathway_name_genes <- c(pathway_name_genes, row$gene_name)
      pathway_name_fc <- c(pathway_name_fc, row$log2FoldChange)
    }
  }
  # Ordering by absolute log fc change
  names(pathway_name_fc) <- pathway_name_genes
  pathway_name_fc <- abs(pathway_name_fc)
  
  pathway_name_fc_names <- pathway_name_fc[order(unlist(pathway_name_fc), decreasing=TRUE)]
  names_list <- names(pathway_name_fc_names)
  sig_names_list <- c()
  if(sig_only == TRUE){
  for(i in 1:length(names_list)){
    if(!is.na(results_useful[results_useful$gene_name == names_list[i],]$padj) & results_useful[results_useful$gene_name == names_list[i],]$padj < 0.01){
      sig_names_list <- c(sig_names_list, names_list[i])
    }
  }
    names_list <- sig_names_list
  }
  return(names_list)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

  