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
                "ddimer_scaled")

pheno_txt = c('Urea', 'Chloride', 
              'C-reactive protein', 'Interleukin-6', 'Neutrophil count', 
              'Lymphocyte count', 'Neutrophil-to-lymphocyte ratio', 
              'D-dimers')

pheno_ann_df = data.frame(cbind(pheno_codes, pheno_txt))
names(pheno_ann_df) = c('Codes', 'Text')


#------------#
# FUNCTIONS  #
#------------#
standardize_blood_phenotypes = function(data){
  
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
  
  return(data)
  
}

clean_corrs_df = function(corrs_df){
  
  new_corrs_df = data.frame(matrix(nrow=0, ncol=ncol(corrs_df)))
  names(new_corrs_df) = names(corrs_df)
  
  mirnas = unique(corrs_df$miRNA)
  
  for (mirna in mirnas){
    
    tmp_df = corrs_df %>%
      filter(miRNA == mirna)
    
    genes = unique(tmp_df$Gene)
    
    for (gene in genes){
      
      tmp2_df = tmp_df %>%
        filter(Gene == gene)
      
      top_transcript = tmp2_df %>%
        filter(FDR.P == min(tmp2_df$FDR.P)) %>%
        distinct()
      
      if (nrow(top_transcript) != 0){
        top_transcript = top_transcript[1,]
      }
      
      new_corrs_df = rbind(new_corrs_df, top_transcript)
    }
    
  }
  
  return(new_corrs_df)
}

convert_string_to_list = function(pheno_txt){
  
  pheno_txt_list = trimws(unlist(strsplit(pheno_txt, ",")))
  pheno_codes_list = c()
  
  for (pheno_txt in pheno_txt_list){
    pheno_code = pheno_ann_df[which(pheno_ann_df$Text == pheno_txt), 'Codes']
    pheno_codes_list = c(pheno_codes_list, pheno_code)
  }
  
  return(pheno_codes_list)
}

mediation_analysis = function(mirna, transcript, pheno){
  
  tmp_df = data %>%
    dplyr::select(all_of(mirna), all_of(transcript), all_of(pheno))
  
  names(tmp_df) = c('MIRNA', 'RNA', 'PHENO')
  
  tmp_df = tmp_df %>%
    filter(is.na(MIRNA) == FALSE) %>%
    filter(is.na(RNA) == FALSE) %>%
    filter(is.na(PHENO) == FALSE)
  
  fit_mediator = lm(RNA ~ MIRNA, tmp_df)
  fit_pheno = lm(PHENO ~ RNA + MIRNA, tmp_df)
  
  results = summary(mediate(fit_mediator, fit_pheno, 
                            treat='MIRNA', mediator='RNA', boot=TRUE))
  
  estimate_ACME = round(results$d0, 3)
  ci_ACME = paste(round(results$d0.ci, 3), collapse='-')
  p_ACME = results$d0.p
  
  results_vector = c(mirna, transcript, pheno, estimate_ACME, ci_ACME, p_ACME)
  
  return(results_vector)
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
data = read.csv('~/Desktop/miRNA/rnaseq/all_data_mirna_rna_oasis.csv')

# Load correlated pairs
corrs_df = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_7.csv')
corrs_df = clean_corrs_df(corrs_df)

# Load gene summaries (manually extracted from GeneCard/Google)
gene_summaries = read.csv('~/Desktop/miRNA/oasis/gene_summaries.csv')
gene_summaries[gene_summaries == ''] = NA


#-------------------------------#
# STANDARDIZE BLOOD PHENOTYPES  #
#-------------------------------#
data = standardize_blood_phenotypes(data)


#---------------------------#
# MEDIATION ANALYSIS 2:     #
# miRNA - RNA - phenotype   #
#---------------------------#

## Negative miRNA-RNA correlations
mediation_results_corrs_df = data.frame(matrix(nrow=0, ncol=7))
names(mediation_results_corrs_df) = c('MIRNA', 'RNA', 'PHENO', 
                                      'Mediated effect', 'CI', 'p', 'GENE')

for (i in 1:nrow(corrs_df)){

  mirna = corrs_df[i,'miRNA']
  mirna = str_replace_all(mirna, '-', '.')
  transcript = corrs_df[i,'Transcript.ID']
  gene = trimws(corrs_df[i,'Gene'])
  pheno_text = corrs_df[i,'Associated_phenotypes']
  
  pheno_list = convert_string_to_list(pheno_text)
  
  for (pheno in pheno_list){
    mediation_res = mediation_analysis(mirna, transcript, pheno)
    mediation_results_corrs_df[nrow(mediation_results_corrs_df)+1, ] = 
      c(mediation_res, gene)
  }
  
}

other = mediation_results_corrs_df


# Add BONF correction
bonf_p = 0.05/nrow(mediation_results_corrs_df)
mediation_results_corrs_df = add_sig_stars(mediation_results_corrs_df, 'p', 'Significance', bonf_p)

# Add IPA target confidence
mediation_results_corrs_df$MIRNA = str_replace_all(mediation_results_corrs_df$MIRNA, "\\.", "-")
mediation_results_corrs_df = merge(mediation_results_corrs_df, 
                                   corrs_df %>% dplyr::select(miRNA, Transcript.ID, Confidence),
                                   by.x = c('MIRNA', 'RNA'), 
                                   by.y = c('miRNA', 'Transcript.ID'))


# Add genes 
ordered_df = mediation_results_corrs_df[order(mediation_results_corrs_df$p),]
ordered_df = replace_names(ordered_df, 'PHENO', pheno_codes, pheno_txt)


#-----------------------------------------------#
# SUPPLEMENTARY TABLE 9:                        #
# Mediation analysis: miRNA - RNA - phenotype   #
#-----------------------------------------------#
ordered_df = ordered_df %>%
  dplyr::select(MIRNA, RNA, GENE, PHENO, `Mediated effect`, CI, p, Significance, Confidence) %>%
  filter(Significance == '*')

ordered_df = ordered_df[order(ordered_df$Confidence), ]

for (i in 1:nrow(ordered_df)){
  
  gene = ordered_df[i,'GENE']
  desc = unique(trimws(gene_summaries[gene_summaries$Gene == gene &
                                        is.na(gene_summaries$Gene.summary) == FALSE, 'Gene.summary']))
  
  if (length(desc) != 0){
    ordered_df[i, 'Gene summary'] = desc
  } else {
    print(gene)
  }
  
}

write.csv(ordered_df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_9.csv',
          row.names=FALSE)


## Process mediation table to inform ms text
unique_exp_obs = ordered_df %>%
  filter(Confidence == 'Experimentally Observed') %>%
  dplyr::select(MIRNA, GENE, PHENO) %>%
  distinct()


n_unique_triplets_rna = nrow(ordered_df %>%
                               dplyr::select(MIRNA, RNA, PHENO) %>% distinct())

n_unique_triplets_gene = nrow(ordered_df %>%
                               dplyr::select(MIRNA, GENE, PHENO) %>% distinct())

n_robust_p = nrow(ordered_df %>%
                    filter(Significance == '*'))

n_robust_p_and_exp_obs = nrow(ordered_df %>%
                                filter(Significance == '*') %>%
                                filter(Confidence == 'Experimentally Observed'))

exp_obs_robust_mirna_gene_pheno = ordered_df %>%
  filter(Significance == '*') %>%
  filter(Confidence == 'Experimentally Observed') %>%
  dplyr::select(MIRNA, GENE, PHENO) %>% distinct()

n_unique_mirna_gene_pheno_robust_p_and_exp_obs = nrow(exp_obs_robust_mirna_gene_pheno)

print(str_interp("There are a total of ${n_unique_triplets_rna} unique miRNA-RNA-pheno triplets."))
print(str_interp("There are a total of ${n_unique_triplets_gene} unique miRNA-gene-pheno triplets."))
print(str_interp("There are a total of ${n_robust_p} triplets with robust P (P < ${bonf_p})."))
print(str_interp("There are a total of ${n_unique_mirna_gene_pheno_robust_p_and_exp_obs} triplets with robust P that are exp. observed."))


