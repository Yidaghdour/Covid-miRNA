rm(list=ls())

library(dplyr)
library(mediation)
library(stringr)

#---------#
# LISTS   #
#---------#

pheno_codes = c("urea_scaled", "chloride_scaled",
                "crp_scaled", "il6_scaled", "neutrophil_scaled", 
                "lymphocyte_scaled", "neutro_lympho_ratio_scaled",
                "ddimer_scaled", 'out_ICU')

pheno_txt_og = c('Urea', 'Chloride', 
              'C-reactive protein', 'Interleukin-6', 'Neutrophil count', 
              'Lymphocyte count', 'Neutrophil-to-lymphocyte ratio', 
              'D-dimers', 'ICU admission')

pheno_ann_df = data.frame(cbind(pheno_codes, pheno_txt_og))
names(pheno_ann_df) = c('Codes', 'Text')

#------------#
# FUNCTIONS  #
#------------#

convert_string_to_list = function(pheno_txt){
  
  pheno_txt_list = trimws(unlist(strsplit(pheno_txt, ",")))
  pheno_codes_list = c()
  
  for (pheno_txt in pheno_txt_list){
    pheno_code = pheno_ann_df[which(pheno_ann_df$Text == pheno_txt), 'Codes']
    pheno_codes_list = c(pheno_codes_list, pheno_code)
  }
  
  return(pheno_codes_list)
}

mediation_analysis = function(snp, mirna, pheno){
  
  tmp_df = data %>%
    dplyr::select(all_of(snp), all_of(mirna), all_of(pheno))
  
  names(tmp_df) = c('SNP', 'MIRNA', 'PHENO')
  
  tmp_df = tmp_df %>%
    filter(is.na(SNP) == FALSE) %>%
    filter(is.na(MIRNA) == FALSE) %>%
    filter(is.na(PHENO) == FALSE)
  
  fit_mediator = lm(MIRNA ~ SNP, tmp_df)
  fit_pheno = lm(PHENO ~ MIRNA + SNP, tmp_df)
  
  results = summary(mediate(fit_mediator, fit_pheno, 
                            treat='SNP', mediator='MIRNA', boot=TRUE))
  
  estimate_ACME = round(results$d0, 3)
  ci_ACME = paste(round(results$d0.ci, 3), collapse='-')
  p_ACME = results$d0.p
  
  results_vector = c(snp, mirna, pheno, estimate_ACME, ci_ACME, p_ACME)
  
  return(results_vector)
}

add_fdr_correction = function(df){
  
  df$p = as.numeric(df$p)
  
  df$FDR_P = p.adjust(df$p, 
                      method='fdr', n=nrow(df))
  
  n_sig_nominal = nrow(df %>% filter(p < 0.05))
  n_sig_fdr = nrow(df %>% filter(FDR_P < 0.05))
  
  print(str_interp("There are ${n_sig_nominal} associations with nominal P < 0.05."))
  print(str_interp("There are ${n_sig_fdr} associations with FDR P < 0.05."))
  
  return(df)
  
}

add_sig_stars = function(table_df, p_col_name, significance_col_name, bonf_p){
  
  table_df[,p_col_name] = as.numeric(table_df[,p_col_name])
  
  for (i in 1:nrow(table_df)){
    if (table_df[i, p_col_name] < bonf_p){
      table_df[i, significance_col_name] = '*'
    } else {
      table_df[i, significance_col_name] = 'ns'
    }
  }
  
  return(table_df)
  
}

replace_names = function(df, col_name, values_to_replace, values_to_replace_with){
  
  ann_df = data.frame(cbind(values_to_replace, values_to_replace_with))
  names(ann_df) = c('col_1' ,'col_2')
  
  df[,col_name] = as.character(df[,col_name])
  
  for (i in 1:nrow(df)){
    value = df[i, col_name]
    
    if (value %in% ann_df$col_1){
      replacement_value = ann_df[ann_df$col_1 == value, 'col_2']
      df[i, col_name] = replacement_value
    }
    
  }
  
  return(df)
  
}


#-------------#
# LOAD DATA   #
#-------------#

data = read.csv('~/Desktop/miRNA/mirna_processed/n96/mirna_oasis_5_in_50_percent_standard_clinical.csv')
mirna_geno_data = read.csv('~/Desktop/miRNA/data/relevant_geno_mirna_oasis.csv')

eqtls = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_10.csv')
top_eqtls = read.csv('~/Desktop/covid_paper_outputs_oasis/Table_2.csv')
  

#-------------------------------#
# STANDARDIZE BLOOD PHENOTYPES  #
#-------------------------------#

data$FIRST_UREA_LVL[data$FIRST_UREA_LVL == 0] = NA
data$urea_scaled = scale(data$FIRST_UREA_LVL, center=TRUE, scale=TRUE)

data$FIRST_CHLORIDE_LVL[data$FIRST_CHLORIDE_LVL == 0] = NA
data$chloride_scaled = scale(data$FIRST_CHLORIDE_LVL, center=TRUE, scale=TRUE)

data$FIRST_C_REACTIVE_PROT[data$FIRST_C_REACTIVE_PROT == 0] = NA
data$crp_scaled = scale(data$FIRST_C_REACTIVE_PROT, center=TRUE, scale=TRUE)

data$FIRST_INTERLEUKIN6[data$FIRST_INTERLEUKIN6 == 0] = NA
data$il6_scaled = scale(data$FIRST_INTERLEUKIN6, center=TRUE, scale=TRUE)

data$FIRST_ABSOLUTE_NEUTROPHIL[data$FIRST_ABSOLUTE_NEUTROPHIL == 0] = NA
data$neutrophil_scaled = scale(data$FIRST_ABSOLUTE_NEUTROPHIL, center=TRUE, scale=TRUE)

data$FIRST_ABSOLUTE_LYMPHOCYTE[data$FIRST_ABSOLUTE_LYMPHOCYTE == 0] = NA
data$lymphocyte_scaled = scale(data$FIRST_ABSOLUTE_LYMPHOCYTE, center=TRUE, scale=TRUE)

data$FIRST_NEUTRO_TO_LYMPHO_RATIO[data$FIRST_NEUTRO_TO_LYMPHO_RATIO == 0] = NA
data$neutro_lympho_ratio_scaled = scale(data$FIRST_NEUTRO_TO_LYMPHO_RATIO, center=TRUE, scale=TRUE)

data$FIRST_D_DIMER[data$FIRST_D_DIMER == 0] = NA
data$ddimer_scaled = scale(data$FIRST_D_DIMER, center=TRUE, scale=TRUE)

scaled_vars = c('urea_scaled', 'chloride_scaled', 'crp_scaled', 
                'il6_scaled', 'neutrophil_scaled', 'lymphocyte_scaled', 
                'neutro_lympho_ratio_scaled', 'ddimer_scaled', 'out_ICU')

data_relevant = data %>%
  dplyr::select(Barcode, Age_numeric, Gender_binary, all_of(scaled_vars))

data = merge(data_relevant, mirna_geno_data, by='Barcode')


#-----------------------#
# MEDIATION ANALYSES    #
#-----------------------#

## TOP eQTLS ------------------------------------------------------------
mediation_results_top_eqtls_df = data.frame(matrix(nrow=0, ncol=6))
names(mediation_results_top_eqtls_df) = c('SNP', 'MIRNA', 'PHENO', 'Mediated effect', 'CI', 'p')

for (i in 1:nrow(top_eqtls)){

  snp = top_eqtls[i,'SNP']
  mirna = top_eqtls[i,'MIRNA']
  pheno_txt = paste(top_eqtls[i,'Negative_associations'], 
                    top_eqtls[i,'Positive_associations'], ",")
  
  pheno_list = convert_string_to_list(pheno_txt)
  mirna = str_replace_all(mirna, '-', '.')
  
  for (pheno in pheno_list){
    
    mediation_res = mediation_analysis(snp, mirna, pheno)
    mediation_results_top_eqtls_df[nrow(mediation_results_top_eqtls_df)+1, ] = mediation_res
    
  }
  
}

raw = mediation_results_top_eqtls_df

# Add FDR correction
bonf_p = 0.05/nrow(mediation_results_top_eqtls_df)
mediation_results_top_eqtls_df = add_sig_stars(mediation_results_top_eqtls_df, 'p', 'Significance', bonf_p)

ordered_df = mediation_results_top_eqtls_df[order(mediation_results_top_eqtls_df$p),]
ordered_df = replace_names(ordered_df, 'PHENO', pheno_codes, pheno_txt_og)

ordered_df$MIRNA = str_replace_all(ordered_df$MIRNA, "\\.", "-")

#--------------------------------------------#
# SUPPLEMENTARY TABLE 12:                    #
# Mediations: SNP--miRNA--blood phenotype    #
#--------------------------------------------#
write.csv(ordered_df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_12.csv',
          row.names=FALSE)

