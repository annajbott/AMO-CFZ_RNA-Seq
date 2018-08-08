library(tidyverse)
library("DESeq2")
library("IHW")

##################################
## Importing kallisto gene data ##
##################################

gene_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/genes.tsv.gz")
gene_tib2 <- as.data.frame(gene_tib)[-1]
rownames(gene_tib2) <- as.data.frame(gene_tib)[,1]
head(gene_tib2)

#transcripts_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/transcripts.tsv.gz")

# Looking at data structure
summary(gene_tib2)

# Remove sample number at the end
colnames(gene_tib2) <- sub("_[^_]+$", "", colnames(gene_tib2))
# Remove index number
colnames(gene_tib2) <-gsub("^i[0-9][0-9]-","",colnames(gene_tib2))


# Genes going down, sample name going along
CFZ_genes <- select(gene_tib2, -contains("WT"))

# Separate by time points, might need to move later
CFZ_6hr <- select(CFZ_genes, -contains("t24"))
CFZ_24hr <- select(CFZ_genes, -contains("t6"))

coldata <- read.csv("/home/osboxes/Documents/analysis/AMO-CFZ_RNA-Seq/sample_sheet2.csv", row.names =1)
coldata <- coldata[-1]
names(coldata) <- c("Compound", "TimePoint")

# Replace underscores for hyphens
rownames(coldata) <- gsub("_","-",rownames(coldata))

# Ensure the rows of sample sheet match the columns of the gene matrix
rownames(coldata)
colnames(CFZ_genes)
all(rownames(coldata) %in% colnames(CFZ_genes)) # True 

# They have the same values but not the same order
all(rownames(coldata) == colnames(CFZ_genes)) # False
# Put them in the same order as the sample sheet
CFZ_genes <- CFZ_genes[, rownames(coldata)]
# Test they are in the same order now
all(rownames(coldata) == colnames(CFZ_genes)) # Now true

## Separate into 6 and 24 hour time points ##
## --------------------------------------- ##

CFZ_6hr <- select(CFZ_genes, -contains("t24"))
CFZ_24hr <- select(CFZ_genes, -contains("t6"))

## 6 and 24 hour versions of smaple sheet, keep row names
coldata_6<-coldata[coldata$TimePoint == 6, ]
coldata_24 <- coldata[coldata$TimePoint == 24, ]

# Check 6 and 24 hour filtered tables also match
all(rownames(coldata_6) == colnames(CFZ_6hr)) # TRUE
all(rownames(coldata_24) == colnames(CFZ_24hr)) # TRUE

## Write CSV files of data ##
## ----------------------- ##
write.csv(CFZ_genes, file = "CFZ_genes.csv")
write.csv(CFZ_6hr, file = "CFZ_genes_6hr.csv")
write.csv(CFZ_24hr, file = "CFZ_genes_24hr.csv")

write.csv(coldata, file = "coldata.csv")
write.csv(coldata_6, file = "coldata_6hr.csv")
write.csv(coldata_24, file = "coldata_24hr.csv")


