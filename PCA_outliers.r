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

## Sample distances- heat map ##
## -------------------------- ##
sampleDists <- dist(t(assay(rld24)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld24$Compound
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

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

res_26_6_useful <- res_useful(res_26_lfc, tf = TRUE)
res_22_6_useful <- res_useful(res_22_lfc, tf = TRUE)
res_13_6_useful <- res_useful(res_13_lfc, tf = TRUE)
res_18_6_useful <- res_useful(res_18_lfc, tf = TRUE)
res_26_24_useful <- res_useful(res_26_lfc_24, tf = TRUE)
res_22_24_useful <- res_useful(res_22_lfc_24, tf = TRUE)
res_13_24_useful <- res_useful(res_13_lfc_24, tf = TRUE)
res_18_24_useful <- res_useful(res_18_lfc_24, tf = TRUE)

save(res_26_6_useful, res_22_6_useful, res_13_6_useful, res_18_6_useful,res_26_24_useful ,res_22_24_useful , res_13_24_useful , res_18_24_useful, file = "useful_deseq_results.RData")



######################
## Pathway Analysis ##
######################

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

gs <- unique(unlist(kegg.sigmet[rownames(kegg_26$greater)[1]]))
fc_kegg_26 <- cbind(fc_kegg_26) # From biostar forum
essData <- essGene(gs, fc_kegg_26, ref =NULL, samp =NULL)
convertIDs(rownames(as.data.frame(essData)),"ENTREZID","SYMBOL", org.Hs.eg.db)

# Gene Ontology #
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
kegg_26_6_go <- pathway_full(dds = dds_deseq, contrast = c("Compound", "26", "DMSO"), gset = gobpsets,  same_direction = TRUE)

lapply(kegg_26_6_go[[1]], head)


## XGR ##
## --- ##

library("XGR")
load("useful_deseq_results.RData")
ontology <- "MsigdbC2KEGG"
# Use XGR and deseq2 results for differentially expressed genes in xEnricherGenes function to find key pathways affected


eTerm_26_6 <- enricher_analysis(dds_deseq, c("Compound", "26", "DMSO"), ontology = ontology, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)
eTerm_13_6 <- enricher_analysis(dds_deseq, c("Compound", "13", "DMSO"), result_lfc = res_13_6_useful, ontology = ontology, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)

# Use log fold change threshold rather than q values, as only 2 replicates for DMSO, so SE is not very useful
eTerm_26_24 <- enricher_analysis(dds_deseq24, c("Compound", "26", "DMSO"), ontology = ontology, alpha = 0.05, foldchange_threshold = 1, number_top_genes = 200)
eTerm_13_24 <- enricher_analysis(dds_deseq24, c("Compound", "13", "DMSO"), ontology = ontology, alpha = 0.05, foldchange_threshold = 1, number_top_genes = 200)

eTerm_18_6 <- enricher_analysis(dds_deseq, c("Compound", "18", "DMSO"), ontology = ontology, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)

eTerm_18_6_max <- enricher_analysis(dds_deseq, c("Compound", "18", "DMSO"), ontology = "MsigdbC2REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 2000)

xEnrichViewer(eTerm_26_6)
bp <- xEnrichBarplot(eTerm_26_24, top_num="auto", displayBy="adjp")
print(bp)

list_eTerm <- list(eTerm_26_6, eTerm_26_24, eTerm_13_6, eTerm_13_24)
names(list_eTerm) <- c('NCP26- 6hr', 'NCP26- 24hr', 'MAZ13- 6hr', 'MAZ13- 24hr')
bp_Pathway <- xEnrichCompare(list_eTerm, displayBy="fc", FDR.cutoff=5e-3, wrap.width=50)
bp_Pathway + theme(axis.text.y=element_text(size=10))


eTerm_26_6_re <- enricher_analysis(dds_deseq, c("Compound", "26", "DMSO"), result_lfc = res_26_6_useful, ontology = "MsigdbC2REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200, same.dir =  TRUE)
eTerm_13_6_re <- enricher_analysis(dds_deseq, c("Compound", "13", "DMSO"), result_lfc = res_13_6_useful, ontology = "MsigdbC2REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200, same.dir =  TRUE)
eTerm_13_6_re2 <- enricher_analysis(dds_deseq, c("Compound", "13", "DMSO"), result_lfc = res_13_6_useful, ontology = "REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200, same.dir =  TRUE)

bp <- xEnrichBarplot(eTerm_13_24_re, top_num=20, displayBy="adjp")
print(bp)

# Use log fold change threshold rather than q values, as only 2 replicates for DMSO, so SE is not very useful
eTerm_26_24_re <- enricher_analysis(dds_deseq24, c("Compound", "26", "DMSO"), result_lfc = res_26_24_useful, ontology = "MsigdbC2REACTOME", foldchange_threshold = 0.5,  number_top_genes = 200, same.dir =  TRUE)
eTerm_13_24_re <- enricher_analysis(dds_deseq24, c("Compound", "13", "DMSO"), result_lfc = res_13_24_useful, ontology = "MsigdbC2REACTOME", foldchange_threshold = 0.5, number_top_genes = 200, same.dir =  TRUE)

ATF4_genes <- enriched_genes_annotation(eterm = eTerm_26_6_re, pathway_entrez_array = eTerm_26_6_re$annotation$REACTOME_ACTIVATION_OF_GENES_BY_ATF4, results_useful = res_26_6_useful, sig_only = TRUE)
paste(ATF4_genes, collapse = ", ")
UPR_genes <- enriched_genes_annotation(eterm = eTerm_26_6_re, pathway_entrez_array = eTerm_26_6_re$annotation$REACTOME_UNFOLDED_PROTEIN_RESPONSE, results_useful = res_26_6_useful, sig_only = TRUE)
paste(UPR_genes, collapse = ", ")

res_26_6_useful[res_26_6_useful$gene_name == ATF4_genes[5],]

filter(res_26_6_useful, str_detect(gene_name, "SERP1"))

list_eTerm_re <- list(eTerm_26_6_re, eTerm_26_24_re, eTerm_13_6_re, eTerm_13_24_re)
names(list_eTerm_re) <- c('NCP26- 6hr', 'NCP26- 24hr', 'MAZ13- 6hr', 'MAZ13- 24hr')
bp_Pathway_re <- xEnrichCompare(list_eTerm_re, displayBy="fc", FDR.cutoff=5e-4, wrap.width=50)
bp_Pathway_re + theme(axis.text.y=element_text(size=10))

load("useful_deseq_results.RData")

subnet_26_6 <- subneter_analysis(dds_deseq, c("Compound", "26", "DMSO"), result_lfc = res_26_6_useful, alpha = 0.05, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_26_6, pattern=-log10(as.numeric(V(subnet_26_6)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
subnet_26_24 <- subneter_analysis(dds_deseq24, c("Compound", "26", "DMSO"), result_lfc = res_26_24_useful, alpha = 1, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_26_24, pattern=-log10(as.numeric(V(subnet_26_24)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
subnet_13_6 <- subneter_analysis(dds_deseq, c("Compound", "13", "DMSO"), result_lfc = res_13_6_useful, alpha = 0.05, number_top_genes = 500,  subnet.size = 50)
xVisNet(g=subnet_13_6, pattern=-log10(as.numeric(V(subnet_13_6)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
subnet_13_24 <- subneter_analysis(dds_deseq24, c("Compound", "13", "DMSO"), result_lfc = res_13_24_useful, alpha = 1, number_top_genes = 500, subnet.size = 50)
xVisNet(g=subnet_13_24, pattern=-log10(as.numeric(V(subnet_13_24)$significance)), colorbar = FALSE, vertex.shape="circle", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.7, vertex.label.cex = 0.9, vertex.label.color = "black")
#
