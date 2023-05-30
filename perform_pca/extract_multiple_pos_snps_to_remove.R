rm(list=ls())

library(dplyr)
library(stringr)

log = read.table("/scratch/tg1407/COVID_miRNA/data/genotyping_processed/mirna_and_1000g.log", header = FALSE, sep = "\n", fill=TRUE)
names(log) = 'text'

log = log %>%
  filter(str_detect(text, "Warning: Multiple"))

for (i in 1:nrow(log)){

  txt = as.character(log[i,'text'])
  snp_dot = unlist(strsplit(txt, " +"))[7]
  snp = substr(snp_dot, 1, nchar(snp_dot)-1)

  log[i,'SNP'] = snp
}

## Load the other to be exluded SNPs file
snps_to_exclude = read.delim2('/scratch/tg1407/COVID_miRNA/data/genotyping_processed/mirna_and_1000g-merge.missnp', header=FALSE)
names(snps_to_exclude) = 'SNP'

snps_to_exclude_df = data.frame(c(as.character(log$SNP), as.character(snps_to_exclude$SNP)))

snps_to_exclude_df = snps_to_exclude_df %>% distinct()

write.table(snps_to_exclude_df, '/scratch/tg1407/COVID_miRNA/data/genotyping_processed/mirna_and_1000g_snps_to_remove.txt', col.names=FALSE,
            quote=FALSE, row.names=FALSE)
