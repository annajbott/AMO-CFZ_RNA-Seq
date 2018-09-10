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

# Load CFZ genes
CFZ_genes <- read.csv("CFZ_genes.csv", row.names = 1, check.names = FALSE)
coldata_cfz <- read.csv("coldata.csv", row.names = 1)

# Load WT
WT_genes <- read.csv("WT_genes.csv", row.names = 1, check.names = FALSE)
coldata_wt <- read.csv("coldata_wt.csv", row.names = 1)

# Combine data
coldata_DMSO <- rbind(coldata_cfz[coldata_cfz$Compound == "DMSO",],coldata_wt[coldata_wt$Compound == "DMSO",])
coldata_DMSO$CellType <- rep(c("CFZ","WT"), each = 6, len = 12)
coldata_DMSO$CellType <- as.factor(coldata_DMSO$CellType)
coldata_DMSO$TimePoint <- as.factor(coldata_DMSO$TimePoint)
coldata_DMSO <- dplyr::select(coldata_DMSO, -contains("Compound"))

genes_DMSO <- cbind(dplyr::select(CFZ_genes, contains("DMSO")), dplyr::select(WT_genes, contains("DMSO")))
# Check in same order
all(colnames(genes_DMSO) == rownames(coldata_DMSO)) # True, same order


## Using DESeq2 Package ##
## -------------------- ##

# Use cell and timepoint in design

dds_all <-  DESeqDataSetFromMatrix(countData = genes_DMSO,
                               colData = coldata_DMSO,
                               design = CellType~TimePoint)
# Prefiltering, removing rows with sum of counts under 10 #
keep <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[keep,]
dds_all$CellType <- relevel(dds_all$CellType, ref = "WT")

dds_all <- DESeq(dds_all)

# PCA 
vst_all <- vst(dds_all, blind=FALSE)
rld_all <- rlog(dds_all, blind=FALSE)
plotPCA(vst_all, intgroup=c("CellType","TimePoint"))
plotPCA(rld_all, intgroup=c("CellType","TimePoint"))

# heat map
sampleDists <- dist(t(assay(vst_all)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vst_all)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## Separate time points ##
## -------------------- ##
coldata_DMSO
head(genes_DMSO)

coldata_DMSO_6 <- coldata_DMSO[coldata_DMSO$TimePoint == 6,]
coldata_DMSO_24 <- coldata_DMSO[coldata_DMSO$TimePoint == 24,]

genes_DMSO_6 <- dplyr::select(genes_DMSO, -contains("t24"))
genes_DMSO_24 <- dplyr::select(genes_DMSO, -contains("t6"))

## Deseq2 with separate time points ##
dds_dmso_6 <- dds_get(gene_data = genes_DMSO_6, coldata = coldata_DMSO_6, prefilter = 10, reference = "WT", cell_type = TRUE)
dds_dmso_24 <- dds_get(gene_data = genes_DMSO_24, coldata = coldata_DMSO_24, prefilter = 10, reference = "WT", cell_type = TRUE)

res_dsmo_6 <- results_process(dds_dmso_6, contrast = c("CellType", "CFZ", "WT"), alpha = 0.05, significant_only = FALSE)
res_dsmo_24 <- results_process(dds_dmso_24, contrast = c("CellType", "CFZ", "WT"), alpha = 0.05, significant_only = FALSE)

eTerm_dmso_6 <- enricher_analysis(dds_dmso_6, c("CellType", "CFZ", "WT"), ontology = "MsigdbC2REACTOME", result_lfc = res_dsmo_6, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)

eTerm_dmso_24 <- enricher_analysis(dds_dmso_24, c("CellType", "CFZ", "WT"), ontology = "MsigdbC2REACTOME", result_lfc = res_dsmo_24, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)

bp_6 <- xEnrichBarplot(eTerm_dmso_6, top_num="auto", displayBy="adjp")
print(bp_6)
# Bar plot for 24 hour
bp_24 <- xEnrichBarplot(eTerm_dmso_24, top_num="auto", displayBy="adjp")
print(bp_24)

list_eTerm <- list(eTerm_dmso_6, eTerm_dmso_24)
names(list_eTerm) <- c('CFZ vs WT (DMSO)- 6hr', 'CFZ vs WT (DMSO)- 24hr' )
bp_Pathway <- xEnrichCompare(list_eTerm, displayBy="adjp", FDR.cutoff=5e-2, wrap.width=50)
bp_Pathway + theme(axis.text.y=element_text(size=10))



subnet_dmso_6 <- subneter_analysis(dds_dmso_6, c("CellType", "CFZ", "WT"), ontology = "MsigdbC2REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500)
xVisNet(g=subnet_dmso_6, pattern=-log10(as.numeric(V(subnet_dmso_6)$significance)),vertex.shape="sphere", colormap="yr", signature = FALSE, newpage = FALSE)
#

# Top genes and TFs
# TFs
TF_csv <- read_csv("temp/TFCheckpoint.csv")

colnames(TF_csv)
hs_tf <- dplyr::select(TF_csv, c("gene_symbol", "entrez_human", "gene_name", "synonym", "DbTF"))
head(hs_tf)
#sig genes

sig_dmso_6 <- res_ouput(res_dsmo_6, gset = NULL, kegg_output = FALSE, alpha = 0.05)
sig_dmso_6 <- sig_dmso_6[abs(sig_dmso_6$log2FoldChange) >= 1,]
dmso_6_tf <- find_tf(sig_dmso_6$gene_name, hs_tf$gene_symbol)


sig_dmso_24 <- res_ouput(res_dsmo_24, gset = NULL, kegg_output = FALSE, alpha = 0.05)
sig_dmso_24 <- sig_dmso_24[abs(sig_dmso_24$log2FoldChange) >= 1,]
dmso_24_tf <- find_tf(sig_dmso_24$gene_name, hs_tf$gene_symbol)

intersect(dmso_6_tf,dmso_24_tf)


## 6hr vs 24hr ##
## ----------- ##

# Make sure 6 hr and 24 hr match more or less
all_res_dsmo_6 <- res_dsmo_6[!(is.na(res_dsmo_6$padj) | res_dsmo_6$padj==""), ]
all_res_dsmo_24 <- res_dsmo_24[!(is.na(res_dsmo_24$padj) | res_dsmo_24$padj==""), ]

kept <- all_res_dsmo_6[rownames(all_res_dsmo_6) %in% rownames(all_res_dsmo_24),]
all(rownames(all_res_dsmo_24) %in% rownames(kept))
kept2 <- all_res_dsmo_24[rownames(all_res_dsmo_24) %in% rownames(kept),]
# check in the right order
all(rownames(kept) == rownames(kept2)) # Right order

# 6hr vs 24 hr DMSO plot (expect all to be along ab line)
time_vs_dmso <- as.data.frame(cbind(kept$log2FoldChange, kept2$log2FoldChange))
rownames(time_vs_dmso) <- rownames(kept)
colnames(time_vs_dmso) <- c("log_fold_change_6hr", "log_fold_change_24hr")
time_vs_dmso$threshold_status <- ifelse(time_vs_dmso$log_fold_change_6hr >= 1 & time_vs_dmso$log_fold_change_24hr >= 1 | time_vs_dmso$log_fold_change_6hr <= -1 & time_vs_dmso$log_fold_change_24hr <= -1 , TRUE, FALSE)
ggplot(time_vs_dmso, aes(x=log_fold_change_6hr, y=log_fold_change_24hr, color = threshold_status)) + geom_point() +
  labs(x = "Log 2 Fold Change WT vs CFZ (DMSO)- 6hr", y =  "Log 2 Fold Change WT vs CFZ (DMSO)- 24hr", color = "Exceed threshold" ) +
  scale_colour_manual(values = c("Black", "Red")) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


time_vs_dmso$average_fc <- (time_vs_dmso$log_fold_change_6hr + time_vs_dmso$log_fold_change_24hr)/2
###################

# DESeq2
dds_624 <-  DESeqDataSetFromMatrix(countData = genes_DMSO,
                                   colData = coldata_DMSO,
                                   design = ~ TimePoint + CellType)
# Prefiltering, removing rows with sum of counts under 10 #
keep <- rowSums(counts(dds_624)) >= 10
dds_624 <- dds_624[keep,]
dds_624$CellType <- relevel(dds_624$CellType, ref = "WT")

dds_all_624 <- DESeq(dds_624)

# Results
res_dmso_624 <- results_process(dds_all_624, contrast = c("CellType", "CFZ", "WT"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = 0)
res_dmso_624_useful <- res_useful(res_dmso_624, tf= TRUE, alpha = 0.05)

# XGR
eTerm_dmso_ms <- enricher_analysis(dds_all_624, c("CellType", "CFZ", "WT"), ontology = "MsigdbC2REACTOME", result_lfc = res_dmso_624, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)
eTerm_dmso <- enricher_analysis(dds_all_624, c("CellType", "CFZ", "WT"), ontology = "REACTOME", result_lfc = res_dmso_624, alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 200)

# Bar plot
bp_dmso <- xEnrichBarplot(eTerm_dmso, top_num=10, displayBy="adjp", signature =FALSE)
print(bp_dmso)

# Sub network
subnet_dmso <- subneter_analysis(dds_all_624, c("CellType", "CFZ", "WT"), ontology = "MsigdbC2REACTOME", alpha = 0.05, foldchange_threshold = FALSE, number_top_genes = 500, subnet.size = 25)
xVisNet(g=subnet_dmso, pattern=-log10(as.numeric(V(subnet_dmso)$significance)),vertex.shape="sphere", colormap="yr", signature = FALSE, newpage = FALSE, vertex.label.dist =0.4)
#

# compare AMO-CFZ dmso (/ vs WT) with WT treated with CFZ
# Do MA plot

################################
## Genes in enriched pathways ##
################################
# REACTOME Molecular signatures
## G alpa signalling
REACTOME_G_ALPHA_I_SIGNALLING_EVENTS_fc_names <- enriched_genes_annotation(eterm = eTerm_dmso_ms, pathway_entrez_array = eTerm_dmso_ms$annotation$REACTOME_G_ALPHA_I_SIGNALLING_EVENTS, results_useful = res_dmso_624_useful)
paste(REACTOME_G_ALPHA_I_SIGNALLING_EVENTS_fc_names, collapse = ", ")

## TCR signalling
REACTOME_TCR_SIGNALING_fc_names <- enriched_genes_annotation(eterm = eTerm_dmso_ms, pathway_entrez_array = eTerm_dmso_ms$annotation$REACTOME_TCR_SIGNALING, results_useful = res_dmso_624_useful)
paste(REACTOME_TCR_SIGNALING_fc_names, collapse = ", ")

## Interferon gamma signalling 
REACTOME_INTERFERON_GAMMA_SIGNALING_fc_names <- enriched_genes_annotation(eterm = eTerm_dmso_ms, pathway_entrez_array = eTerm_dmso_ms$annotation$REACTOME_INTERFERON_GAMMA_SIGNALING, results_useful = res_dmso_624_useful)
paste(REACTOME_INTERFERON_GAMMA_SIGNALING_fc_names, collapse = ", ")

## MHC calss II antigen presentation
REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION_fc_names <- enriched_genes_annotation(eterm = eTerm_dmso_ms, pathway_entrez_array = eTerm_dmso_ms$annotation$REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION, results_useful = res_dmso_624_useful)
paste(REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION_fc_names, collapse = ", ")

## Semaphorin interactions
REACTOME_SEMAPHORIN_INTERACTIONS_fc_names <- enriched_genes_annotation(eterm = eTerm_dmso_ms, pathway_entrez_array = eTerm_dmso_ms$annotation$REACTOME_SEMAPHORIN_INTERACTION, results_useful = res_dmso_624_useful)
paste(REACTOME_SEMAPHORIN_INTERACTIONS_fc_names, collapse = ", ")

# REACTOME pathways #
## immunoregulatory interactions between a lymphoid and non-lymphoid cell
immuno_lymph_names <- enriched_genes_annotation(eterm = eTerm_dmso, pathway_entrez_array = eTerm_dmso$annotation$`R-HSA-198933`, results_useful = res_dmso_624_useful)
paste(immuno_lymph_names, collapse = ", ")

## Immune system 
immune_system_names <- enriched_genes_annotation(eterm = eTerm_dmso, pathway_entrez_array = eTerm_dmso$annotation$`R-HSA-168256`, results_useful = res_dmso_624_useful)
paste(immune_system_names, collapse = ", ")

## Neutrophil degranulation
neutrophil_degranulation_names <- enriched_genes_annotation(eterm = eTerm_dmso, pathway_entrez_array = eTerm_dmso$annotation$`R-HSA-6798695`, results_useful = res_dmso_624_useful)
paste(neutrophil_degranulation_names, collapse = ", ")

## Cytokine Signaling in Immune system
cytokine_signalling_immune_system_names <- enriched_genes_annotation(eterm = eTerm_dmso, pathway_entrez_array = eTerm_dmso$annotation$`R-HSA-1280215`, results_useful = res_dmso_624_useful)
paste(cytokine_signalling_immune_system_names, collapse = ", ")

## PI3K/AKT Signaling in Cancer
PI3K_AKT_names <- enriched_genes_annotation(eterm = eTerm_dmso, pathway_entrez_array = eTerm_dmso$annotation$`R-HSA-194315`, results_useful = res_dmso_624_useful)
paste(PI3K_AKT_names, collapse = ", ")
