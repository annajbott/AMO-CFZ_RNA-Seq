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

# Looking at data structure
summary(gene_tib2)

# Remove sample number at the end
colnames(gene_tib2) <- sub("_[^_]+$", "", colnames(gene_tib2))
# Remove index number
colnames(gene_tib2) <-gsub("^i[0-9][0-9]-","",colnames(gene_tib2))


# Genes going down, sample name going along
CFZ_genes <- dplyr::select(gene_tib2, -contains("WT"))

# Separate by time points, might need to move later
CFZ_6hr <- dplyr::select(CFZ_genes, -contains("t24"))
CFZ_24hr <- dplyr::select(CFZ_genes, -contains("t6"))

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

CFZ_6hr <- dplyr::select(CFZ_genes, -contains("t24"))
CFZ_24hr <- dplyr::select(CFZ_genes, -contains("t6"))

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

################
## WT samples ##
################

# Genes going down, sample name going along
WT_genes <- dplyr::select(gene_tib2, -contains("CFZ"))

# Separate by time points, might need to move later
WT_6hr <- dplyr::select(WT_genes, -contains("t24"))
WT_24hr <- dplyr::select(WT_genes, -contains("t6"))

coldata_wt <- read.csv("/home/osboxes/Documents/analysis/AMO-CFZ_RNA-Seq/sample_sheet_wt.csv", row.names =1)
coldata_wt <- coldata[c("Compound", "TimePoint")]
names(coldata_wt) <- c("Compound", "TimePoint")

# Replace underscores for hyphens
rownames(coldata_wt) <- gsub("_","-",rownames(coldata_wt))

# Ensure the rows of sample sheet match the columns of the gene matrix
all(rownames(coldata_wt) %in% colnames(WT_genes)) # True 

# They have the same values but not the same order
all(rownames(coldata_wt) == colnames(WT_genes)) # True but lets make sure anyway
# Put them in the same order as the sample sheet
WT_genes <- WT_genes[, rownames(coldata_wt)]
# Test they are in the same order now
all(rownames(coldata_wt) == colnames(WT_genes)) # Now true

## Separate into 6 and 24 hour time points ##
## --------------------------------------- ##

# Separate by time points, might need to move later
WT_6hr <- dplyr::select(WT_genes, -contains("t24"))
WT_24hr <- dplyr::select(WT_genes, -contains("t6"))

## 6 and 24 hour versions of smaple sheet, keep row names
coldata_wt_6 <-coldata_wt[coldata_wt$TimePoint == 6, ]
coldata_wt_24 <- coldata_wt[coldata_wt$TimePoint == 24, ]

# Check 6 and 24 hour filtered tables also match
all(rownames(coldata_wt_6) == colnames(WT_6hr)) # TRUE
all(rownames(coldata_wt_24) == colnames(WT_24hr)) # TRUE

## Write CSV files of data ##
## ----------------------- ##
write.csv(WT_genes, file = "WT_genes.csv")
write.csv(WT_6hr, file = "WT_genes_6hr.csv")
write.csv(WT_24hr, file = "WT_genes_24hr.csv")

write.csv(coldata_wt, file = "coldata_wt.csv")
write.csv(coldata_wt_6, file = "coldata_wt_6hr.csv")
write.csv(coldata_wt_24, file = "coldata_wt_24hr.csv")
