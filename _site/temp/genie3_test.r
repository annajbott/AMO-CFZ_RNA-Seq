library(GENIE3)
library(tidyverse)

head(CFZ_6hr)
keep <- c()
for(i in 1:length(rownames(CFZ_6hr))){
  if(sum(CFZ_6hr[i,]) > 20){
    keep <- c(keep,i)
  }
}
CFZ_6hr_filt <- CFZ_6hr[keep,]
weightMat <- GENIE3(as.matrix(CFZ_6hr_filt))

#

load("useful_deseq_results.RData")


samples_list <- list(res_26_6_useful, res_26_24_useful,res_22_6_useful, res_22_24_useful,res_13_6_useful, res_13_24_useful,res_18_6_useful,res_18_24_useful)
gene_list <- c("EIF4EBP1", "ASNS", "GPT2", "DDIT3" , "ATF4", "ATF3","SLC6A9", "TRIB3")
plot_table <- data.frame()

plot_table[1:2,1] <- "26"
plot_table[3:4,1] <- "22"
plot_table[5:6,1] <- "13"
plot_table[7:8,1] <- "18"
plot_table[c(1,3,5,7),2] <- "6"
plot_table[c(2,4,6,8),2] <- "24"
rownames(plot_table) <- c("26_6", "26_24", "22_6", "22_24", "13_6", "13_24", "18_6", "18_24")
colnames(plot_table) <- c("Compound", "TimePoint")


for(i in seq_along(samples_list)){
  for(j in seq_along(gene_list)){
    subset_data <- samples_list[[i]] %>% filter(str_detect(gene_name, gene_list[j]))
    plot_table[i,j + 2] <- subset_data[1,"log2FoldChange"]
  }
}
colnames(plot_table) <- c("Compound", "TimePoint",gene_list)
plot_table$Compound <- as.factor(plot_table$Compound)
plot_table$TimePoint <- as.factor(plot_table$TimePoint)

plot_table_tidy <- gather(plot_table,"gene", "Log2FoldChange", -c(Compound, TimePoint))
plot_table_tidy$Compound <- factor(plot_table_tidy$Compound,levels = c("13","26","18","22"))

ggplot(data=plot_table_tidy, aes(x=gene, y=Log2FoldChange, fill = Compound)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge()) + 
  labs(x = "Gene", y =  "Log 2 Fold Change")

## AARs genes
gene_list2 <- c("GARS", "WARS", "CARS", "TARS", "AARS", "SARS", "YARS", "IARS", "NARS", "MARS")
plot_table <- data.frame()
plot_table[1:2,1] <- "26"
plot_table[3:4,1] <- "22"
plot_table[5:6,1] <- "13"
plot_table[7:8,1] <- "18"
plot_table[c(1,3,5,7),2] <- "6"
plot_table[c(2,4,6,8),2] <- "24"
rownames(plot_table) <- c("26_6", "26_24", "22_6", "22_24", "13_6", "13_24", "18_6", "18_24")
colnames(plot_table) <- c("Compound", "TimePoint")


for(i in seq_along(samples_list)){
  for(j in seq_along(gene_list2)){
    subset_data <- samples_list[[i]] %>% filter(str_detect(gene_name, gene_list2[j]))
    plot_table[i,j + 2] <- subset_data[1,"log2FoldChange"]
  }
}
colnames(plot_table) <- c("Compound", "TimePoint",gene_list2)
plot_table$Compound <- as.factor(plot_table$Compound)
plot_table$TimePoint <- as.factor(plot_table$TimePoint)

plot_table_tidy <- gather(plot_table,"gene", "Log2FoldChange", -c(Compound, TimePoint))
plot_table_tidy$Compound <- factor(plot_table_tidy$Compound,levels = c("13","26","18","22"))

ggplot(data=plot_table_tidy, aes(x=gene, y=Log2FoldChange, fill = Compound)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge()) + 
  labs(x = "Gene", y =  "Log 2 Fold Change")

