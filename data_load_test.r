library(tidyverse)
library("DESeq2")

gene_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/genes.tsv.gz")

transcripts_tib <- read_tsv("/home/osboxes/Documents/analysis/quant/kallisto.dir/transcripts.tsv.gz")

# Looking at data structure
glimpse(gene_tib)
summary(gene_tib)
for(col_name in colnames(gene_tib)){
  if("WT" in col_name){
    print(col_name)
  }
}
# Genes going down, sample name going along
select