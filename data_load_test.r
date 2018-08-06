library(tidyverse)
library("DESeq2")
library(readxl)

gene_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/genes.tsv.gz")

transcripts_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/transcripts.tsv.gz")

# Looking at data structure
glimpse(gene_tib)
summary(gene_tib)

# Remove sample number at the end
colnames(gene_tib) <- sub("_[^_]+$", "", colnames(gene_tib))
# Remove index number
colnames(gene_tib) <-gsub("^[^-]*-","",colnames(gene_tib))


# Genes going down, sample name going along
CFZ_genes <- select(gene_tib, -contains("WT"))

# Separate by time points
CFZ_6hr <- select(CFZ_genes, -contains("t24"))
CFZ_24hr <- select(CFZ_genes, -contains("t6"))

# Read in sample sheet
coldata <- read_excel("sample_sheet.xlsx")
coldata_6 <- filter(coldata, TimePoint == 6)
coldata_24 <- filter(coldata, TimePoint == 24)




