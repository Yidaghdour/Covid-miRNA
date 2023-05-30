rm(list=ls())

library(dplyr)
library(stringr)

setwd('~/Downloads')

# Genotyping data
geno_ids = read.table('mirna_91_covid_mirna_bia_maf01_hwe005_geno50_mind50.fam')
geno_ids = geno_ids[,1:2]
names(geno_ids) = c('FID', 'IID')

# miRNA data
mirna_data = read.csv('~/Desktop/COVID_miRNA/Processed_data/mirna_encode_5in50_standard_clinical.csv')
mirna_data = mirna_data %>%
  select(Barcode, MIRLET7E:AC091053.2)

# Covariates data
covar_data = read.csv('~/Desktop/COVID_miRNA/Processed_data/mirna_encode_5in50_standard_clinical.csv')
covar_data = covar_data %>%
  select(Barcode, Age_numeric, Gender_binary)

# Annotation
ann = read.csv('annotation_geno_covid_id.csv')


mirna_list = names(mirna_data)[2:ncol(mirna_data)]

#---------------------------#
# MAKE & SAVE PHENO FILES   #
#---------------------------#
setwd('~/Downloads/mirna_pheno')

for (mirna in mirna_list){

  tmp_data = mirna_data %>%
    select(Barcode, all_of(mirna))

  # Annotate miRNA data
  mirna_data_merged = merge(tmp_data, ann, by.x='Barcode', by.y='COVID_ID')

  mirna_data_merged = mirna_data_merged[,c(3, 3, 2)]
  names(mirna_data_merged) = c('FID', 'IID', 'Pheno')

  write.table(mirna_data_merged, str_interp('pheno_${mirna}.txt'),
              sep='\t', quote=FALSE, row.names=FALSE)
}


#---------------------------#
# MAKE AND SAVE COVAR FILE  #
#---------------------------#

covar_data_merged = merge(covar_data, ann, by.x='Barcode', by.y='COVID_ID')
covar_data_merged = covar_data_merged[,c(4, 4, 2, 3)]
names(covar_data_merged) = c('FID', 'IID', 'AGE', 'GENDER')

write.table(covar_data_merged, 'covars_age_gender.txt',
            sep='\t', quote=FALSE, row.names=FALSE)

#---------------------#
# MAKE A MIRNA LIST   #
#---------------------#
write.table(mirna_list, 'mirna_list_encode.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)


