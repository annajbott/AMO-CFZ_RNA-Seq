source("fun.R")
library(tidyverse)
library("DESeq2")
library(org.Hs.eg.db)
library(VennDiagram)

## Venn Diagram ##
## ------------ ##

load("data_objects.RData")

res_26_lfc <- results_process(dds_deseq, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = TRUE, foldchange_threshold = 1 )
res_13_lfc <- results_process(dds_deseq, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = TRUE, foldchange_threshold = 1 )

# Taken out a replicate so padj values are using sd with only two replicates, instead set threshold as 1 log2 fold change (double expression)
res_26_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "26", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = 1)
res_13_lfc_24 <- results_process(dds_deseq24, contrast = c("Compound", "13", "DMSO"), alpha = 0.05, significant_only = FALSE, foldchange_threshold = 1 )

names_26_6 <- convertIDs(rownames(res_26_lfc),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
names_13_6 <- convertIDs(rownames(res_13_lfc),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
names_26_6 <- names_26_6[!is.na(names_26_6)]
names_13_6 <- names_13_6[!is.na(names_13_6)]

names_26_24 <- convertIDs(rownames(res_26_lfc_24),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
names_13_24 <- convertIDs(rownames(res_13_lfc_24),"ENSEMBL", "SYMBOL", org.Hs.eg.db)
names_26_24 <- names_26_24[!is.na(names_26_24)]
names_13_24 <- names_13_24[!is.na(names_13_24)]


# Generate plot
v <- venn.diagram(list(NCP_26_6hr=names_26_6, MAZ_13_6hr=names_13_6),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 0.8, cex=0.5,
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v)
# See which label is which
lapply(v, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in MAZ13 only 
v[[5]]$label  <- " " #paste(setdiff(names_13_6, names_26_6)[1:15], setdiff(names_13_6, names_26_6)[16:30],setdiff(names_13_6, names_26_6)[31:45],setdiff(names_13_6, names_26_6)[46:60], setdiff(names_13_6, names_26_6)[61:70], collapse="\n")  
# in NCP26 only
v[[6]]$label <- " " # paste(setdiff(names_26_6, names_13_6)[1:7], setdiff(names_26_6, names_13_6)[8:14], collapse="\n")  
# intesection
v[[7]]$label <- " " # paste(intersect(names_13_6, names_26_6)[1:25], intersect(names_13_6, names_26_6)[26:50],intersect(names_13_6, names_26_6)[51:75],intersect(names_13_6, names_26_6)[76:100], intersect(names_13_6, names_26_6)[101:106],collapse = "\n")

#plot  
grid.newpage()
grid.draw(v)


# Generate plot
v24 <- venn.diagram(list(NCP_26_24hr=names_26_24, MAZ_13_24hr=names_13_24),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 0, cex=0.5,
                  filename=NULL)

grid.newpage()
grid.draw(v24)
# See which label is which
lapply(v24, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in NCP 26 only
v24[[5]]$label  <- " " 
# in MAZ 13 only
v24[[6]]$label <- " "  
# intesection
v24[[8]]$label <- " " 
#plot  
grid.newpage()
grid.draw(v24)

#################################
## Check transcription factors ##
#################################

TF_csv <- read_csv("temp/TFCheckpoint.csv")

colnames(TF_csv)
hs_tf <- dplyr::select(TF_csv, c("gene_symbol", "entrez_human", "gene_name", "synonym", "DbTF"))
head(hs_tf)

all_genes <- unique(c(names_13_6,names_13_24,names_26_6,names_26_24))
maz_only_6 <- setdiff(names_13_6, names_26_6)
ncp_only_6 <- setdiff(names_26_6, names_13_6)
both_6 <- intersect(names_26_6, names_13_6)
maz_only_24 <- setdiff(names_13_24, names_26_24)
ncp_only_24 <- setdiff(names_26_24, names_13_24)
both_24 <- intersect(names_26_24, names_13_24)

# Use hs_tf$gene_symbol as reference for function
tf_maz_only_6 <- find_tf(maz_only_6, hs_tf$gene_symbol)
tf_maz_only_24 <- find_tf(maz_only_24, hs_tf$gene_symbol) # No TF.
tf_ncp_only_6 <- find_tf(ncp_only_6, hs_tf$gene_symbol)
tf_ncp_only_24 <- find_tf(ncp_only_24, hs_tf$gene_symbol)
tf_both_6 <- find_tf(both_6, hs_tf$gene_symbol)
tf_both_24 <- find_tf(both_24, hs_tf$gene_symbol)
tf_all <- find_tf(all_genes, hs_tf$gene_symbol)

tf_names_13_6 <- find_tf(names_13_6, hs_tf$gene_symbol)
tf_names_13_24 <- find_tf(names_13_24, hs_tf$gene_symbol)
tf_names_26_6 <- find_tf(names_26_6, hs_tf$gene_symbol)
tf_names_26_24 <- find_tf(names_26_24, hs_tf$gene_symbol)

v_tf <- venn.diagram(list(NCP_26_6hr=tf_names_26_6, NCP_26_24hr = tf_names_26_24, MAZ_13_6hr=tf_names_13_6, MAZ_13_24hr = tf_names_13_24),
                  fill = c("orange", "red", "blue", "green"),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 0.8, cex=0.5,
                  filename=NULL)
grid.newpage()
grid.draw(v_tf)
lapply(v_tf, function(i) i$label)
for(i in 9:23){
  v_tf[[i]]$label  <- " " 
}
grid.newpage()
grid.draw(v_tf)
# No experimental backing up TRIB 3 as transcription factor

## Filtering for TFs backed up with experimental evidence for:
#1) regulation of RNA polymerase II 
#2) specific DNA binding activity. 
tf_maz_only_6_evidence <- filter(hs_tf, gene_symbol %in% tf_maz_only_6) %>%
  filter( DbTF == "yes")
tf_ncp_only_24_evidence <- filter(hs_tf, gene_symbol %in% tf_ncp_only_24) %>%
  filter( DbTF == "yes")
tf_both_6_evidence <- filter(hs_tf, gene_symbol %in% tf_both_6) %>%
  filter( DbTF == "yes")
tf_both_24_evidence <- filter(hs_tf, gene_symbol %in% tf_both_24) %>%
  filter( DbTF == "yes")

# May have to use multiple databases to catch all TFs.

# TFCat keeps crashing terrible database
#tfcat_delim <- read_delim("temp/TFCat_delim.txt", delim = " ")
#tfcat_csv <- read_csv("temp/TFCAT.csv", col_types=list("Gene ID" = 'c', "Description" = '_'	, "DBD Classification (Protein Group/Protein Family)" = '_',	"Reviewer" = '_',	"Judgement" = '_',	"Taxa" = '_', "Reviewer Commments" = '_', "PubMed ID" = '_',	"Function" = '_',	"Species"= '_',	"Evidence Strength"= '_',	"PubMed Comments"= '_'))
tfcat_csv <- read.csv("temp/TFCAT.csv")
a <- read_csv('coldata.csv', col_types=list('Compound'='_', TimePoint='c'))

a
# Human tfdb database
tfdb_txt <- read.delim("temp/human_tfdb.txt")

tfdb_maz_13_6 <- find_tf(names_13_6, tfdb_txt$Symbol )
