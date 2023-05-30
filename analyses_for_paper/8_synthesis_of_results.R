rm(list=ls())

library(dplyr)
library(stringr)

setwd('~/Desktop/covid_paper_outputs_oasis/')

# Load data
icu_mirna = read.csv('Supplementary_Table_4.csv')
icu_mirna_rna_cors = read.csv('Supplementary_Table_5.csv')

blood_mirna = read.csv('Supplementary_Table_6.csv')
blood_mirna_rna_cors = read.csv('Supplementary_Table_7.csv')

icu_mirna_p01 = icu_mirna %>% filter(P.value < 0.01)
blood_mirna_p01 = blood_mirna %>% filter(P.value < 0.01)

mediations_mirna_rna_pheno = read.csv('Supplementary_Table_9.csv')
mediations_snp_mirna_pheno = read.csv('Supplementary_Table_11.csv')

top_eqtls = read.csv('~/Desktop/covid_paper_outputs_oasis/Table_2.csv')


# 1. ICU miRNAs (pos and neg)
icu_pos = icu_mirna_p01 %>% filter(Effect > 0)
icu_neg = icu_mirna_p01 %>% filter(Effect < 0)

n_icu_pos = nrow(icu_pos)
n_icu_neg = nrow(icu_neg)

print(str_interp("There are ${n_icu_pos} + associated and ${n_icu_neg} - associatied miRNAs with ICU. "))


# 2. ICU and Blood (shared associations)
n_blood = nrow(blood_mirna_p01)

shared_icu_blood = intersect(icu_mirna_p01$miRNA, blood_mirna_p01$miRNA)
n_shared_icu_blood = length(shared_icu_blood)

print(str_interp("There are ${n_blood} miRNAs assoc. with a blood phenotype."))
print(str_interp("${n_shared_icu_blood} are shared with ICU."))
print(str_interp("Those are ${shared_icu_blood}"))


# 3. 