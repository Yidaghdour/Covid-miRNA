rm(list=ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(data.table)


#---------#
# LISTS   #
#---------#
ipa_results_dir = '~/Desktop/miRNA/ipa_results/oasis/'

# Blood phenotype codes
pheno_codes = c("urea_scaled", "chloride_scaled",
                "crp_scaled", "il6_scaled", "neutrophil_scaled", 
                "lymphocyte_scaled", "neutro_lympho_ratio_scaled",
                "ddimer_scaled")

pheno_txt = c('Urea', 'Chloride', 
              'C-reactive protein', 'Interleukin-6', 'Neutrophil count', 
              'Lymphocyte count', 'Neutrophil-to-lymphocyte ratio', 
              'D-dimers')

#-------------#
# FUNCTIONS   #
#-------------#

# General
add_sig_stars = function(table_df, p_col_name, significance_col_name){
  
  for (i in 1:nrow(table_df)){
    if (table_df[i, p_col_name] < 0.0001){
      table_df[i, significance_col_name] = '****'
    } else if (table_df[i, p_col_name] < 0.001){
      table_df[i, significance_col_name] = '***'
    } else if (table_df[i, p_col_name] < 0.01){
      table_df[i, significance_col_name] = '**'
    } else if (table_df[i, p_col_name] < 0.05){
      table_df[i, significance_col_name] = '*'
    } else if (table_df[i, p_col_name] > 0.05){
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


# IPA analysis
load_ipa_data = function(ipa_results_dir, filename){
  
  setwd(ipa_results_dir)
  
  raw_main_df = readxl::read_xls(str_interp(filename))
  raw_main_df = raw_main_df[-1,]
  
  names(raw_main_df) = c('mirna_id', 'mirna_id_explained', 
                         'source', 'confidence', 'gene', 'pathway')
  
  main_df = raw_main_df %>%
    dplyr::select(mirna_id, confidence, gene, pathway) %>%
    distinct()
  
  # Replace multiple confidence labels with the highest one 
  for (i in 1:nrow(main_df)){
    
    conf = as.character(main_df[i,'confidence'])
    
    if (is.na(conf) == FALSE){
      
      if (grepl('Experimentally Observed', conf, fixed=TRUE)){
        main_df[i,'confidence'] = 'Experimentally Observed'
      } else if (grepl('High (predicted)', conf, fixed=TRUE)){
        main_df[i,'confidence'] = 'High (predicted)'
      } else if (grepl('Moderate (predicted)', conf, fixed=TRUE)){
        main_df[i,'confidence'] = 'Moderate (predicted)'
      }
    }
  }
  
  return(main_df)
}

annotate_mirna_table = function(corrs_table){
  
  corrs_table$predicted_target = NA
  corrs_table$confidence = NA
  corrs_table$gene_pathways = NA
  
  for (i in 1:nrow(corrs_table)){
    
    mirna = trimws(corrs_table[i,'mirna'])
    gene = trimws(corrs_table[i,'ref_gene_name'])
    
    tmp_df = main_df[grepl(mirna, main_df$mirna_id), ]
    
    if (nrow(tmp_df) != 0){
      
      if (gene %in% tmp_df$gene){
        
        confidence = tmp_df[which(tmp_df$gene == gene), 'confidence']
        pathways = tmp_df[which(tmp_df$gene == gene), 'pathway']
        
        corrs_table[i, 'predicted_target'] = 'Yes'
        corrs_table[i, 'confidence'] = confidence 
        corrs_table[i, 'gene_pathways'] = pathways
        
      } else {
        corrs_table[i, 'predicted_target'] = 'No'
      }
      
    } else if (nrow(tmp_df) == 0){
      corrs_table[i, 'predicted_target'] = 'No information'
    }
    
  }
  
  return(corrs_table)
}

add_associated_blood_pheno = function(df, table_biomarkers){
  
  for (i in 1:nrow(df)){
    
    mirna = df[i,'miRNA']
    
    tmp_df = table_biomarkers %>%
      filter(miRNA == mirna) %>%
      select(Dependent.variable)
    
    tmp_df = replace_names(tmp_df, 'Dependent.variable', pheno_codes, pheno_txt)
    
    assoc_phenotypes = paste(tmp_df$Dependent.variable, collapse=', ')
    
    df[i, 'Associated_phenotypes'] = assoc_phenotypes
  }
  
  return(df)
}

process_canonical_pathways = function(path){
  
  pathways = read.delim2(path)
  pathways$X.log.p.value. = as.numeric(pathways$X.log.p.value.)
  pathways$Ratio = as.numeric(pathways$Ratio)
  pathways$P_value = 10^(-pathways$X.log.p.value.)
  
  pathways_ordered_df = pathways[order(pathways$P_value), ]
  pathways_ordered_df = pathways_ordered_df %>%
    filter(P_value < 0.01) %>%
    dplyr::select(Ingenuity.Canonical.Pathways, Ratio, P_value, Molecules)
  
  pathways_ordered_df$Molecules = str_replace_all(pathways_ordered_df$Molecules,
                                                  ",", ", ")
  
  return(pathways_ordered_df)
}

process_diseases_processes = function(path){
  
  # Diseases & functions
  diseases = read.delim2(path)
  diseases$p.value = as.numeric(diseases$p.value)
  
  unique_categories = unique(unlist(str_split(diseases$Categories, ",")))
  
  unique_categories_df = data.frame(matrix(nrow=0, ncol=3))
  names(unique_categories_df) = c('Category', 'P_value', 'N_molecules')
  
  for (cat in unique_categories){
    
    tmp = diseases[diseases$Categories %like% cat,]
    
    p_min = min(tmp$p.value)
    p_max = max(tmp$p.value)
    
    n_max = max(tmp$X..Molecules)
    
    unique_categories_df[nrow(unique_categories_df)+1, ] =
      c(cat, str_interp("${p_min}-${p_max}"), n_max)
  }
  
  unique_categories_df$N_molecules = as.numeric(unique_categories_df$N_molecules)
  
  diseases = c('Cancer', 'Organismal Injury and Abnormalities', 'Infectious Diseases',
               'Renal and Urological Disease', 'Immunological Disease')
  
  functions = c('Cell Death and Survival', 'Cellular Development', 
                'Cellular Growth and Proliferation', 'Cellular Function and Maintenance',
                'Cell Morphology')
  
  diseases_df = unique_categories_df %>% filter(Category %in% diseases)
  diseases_df$Type = 'Diseases and Disorders'
  
  functions_df = unique_categories_df %>% filter(Category %in% functions)
  functions_df$Type = 'Molecular and Cellular Functions'
  
  df = rbind(diseases_df, functions_df)
  
  return(df)
}


#-------------#
# LOAD DATA   #
#-------------#

# Load the correlation results
blood_s_corrs = read.csv('~/Desktop/miRNA/rnaseq/outputs/corrs_spearman_sig_fdr05.csv')

data = read.csv('~/Desktop/miRNA/rnaseq/all_data_mirna_rna_oasis.csv')

table_biomarkers = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_6.csv')
table_biomarkers = table_biomarkers %>%
  filter(P.value < 0.01)


#--------------------------#
# ANNOTATE miRNA TARGETS   #
#--------------------------#

# IPA targets
main_df = load_ipa_data(ipa_results_dir, 'blood_pheno_mirna_targets_updated.xls')

# Annotate miRNA-RNA pairs 
blood_s_corrs$mirna = str_replace_all(blood_s_corrs$mirna,
                                      '\\.', '-')

blood_s_corrs_ann = annotate_mirna_table(blood_s_corrs)

write.csv(blood_s_corrs_ann, '~/Desktop/miRNA/rnaseq/outputs/corrs_spearman_sig_fdr05_bloodpheno_oasis_annotated.csv')


#---------------------------------#
# PROCESS DATA (informs ms text)  #
#---------------------------------#

## NEGATIVELY CORRELATED miRNA-mRNA 
blood_s_corrs_neg_ann = blood_s_corrs_ann %>%
  filter(corr < 0)

blood_s_corrs_neg_ann = add_sig_stars(blood_s_corrs_neg_ann, 'FDR_P', 'Significance')

n_negative_corrs = nrow(blood_s_corrs_neg_ann)
min_neg_cor = round(min(blood_s_corrs_neg_ann$corr), 2)
max_neg_cor = round(max(blood_s_corrs_neg_ann$corr), 2)
mean_neg_cor = round(mean(blood_s_corrs_neg_ann$corr), 2)
sd_neg_cor = round(sd(blood_s_corrs_neg_ann$corr), 2)

print(str_interp(
  "There are a total of ${n_negative_corrs} negatively-correlated miRNA-mRNA pairs."))

print(str_interp(
  "The range was ${min_neg_cor} to ${max_neg_cor}. The mean was ${mean_neg_cor} (SD = ${sd_neg_cor})"))

n_unique_mirna = length(unique(blood_s_corrs_neg_ann$mirna))

print(str_interp("It concerned ${n_unique_mirna} miRNAs."))


## HOW MANY miRNAS HAVE TARGETS FROM IPA? 
n_unique_mirna = length(unique(blood_s_corrs_neg_ann$mirna))

tmp = blood_s_corrs_neg_ann %>%
  filter(predicted_target != 'No information')

n_mirna_w_target_data = length(unique(tmp$mirna))

print(str_interp(
  "We have target data for ${n_mirna_w_target_data} out of ${n_unique_mirna} miRNAs."))


## EXPERIMENTALLY OBSERVED 
exp_observed = blood_s_corrs_neg_ann %>%
  filter(confidence == 'Experimentally Observed')

n_exp_obs = nrow(exp_observed)
n_mirna = length(unique(exp_observed$mirna))
n_genes = length(unique(exp_observed$ref_gene_name))

print(str_interp(
  "There were ${n_exp_obs} experimentally observed pairs, concerning ${n_mirna} miRNAs and ${n_genes} genes."))


## HIGHLY PREDICTED 
highly_predicted = blood_s_corrs_neg_ann %>%
  filter(confidence == 'High (predicted)')

n_exp_obs = nrow(highly_predicted)
n_mirna = length(unique(highly_predicted$mirna))
n_genes = length(unique(highly_predicted$ref_gene_name))

print(str_interp(
  "There were ${n_exp_obs} highly predicted pairs, concerning ${n_mirna} miRNAs and ${n_genes} genes."))


## EXPERIMENTALLY OBSERVED + HIGHLY PREDICTED 
df = blood_s_corrs_neg_ann %>%
  filter(confidence %in% c('Experimentally Observed',
                           'High (predicted)')) %>%
  select(mirna, rna, ref_gene_name, corr, FDR_P, Significance, confidence) %>%
  arrange(confidence, mirna)

names(df) = c('miRNA', 'Transcript ID', 'Gene', 'Correlation coefficient',
              'FDR P', 'Significance', 'Confidence')

df$`Correlation coefficient` = round(df$`Correlation coefficient`, 3)

df = add_associated_blood_pheno(df, table_biomarkers)


#-------------------------------------------------#
# SUPPLEMENTARY TABLE 7:                          #
# Experimentally observed and highly predicted    #
# miRNA-RNA pairs with negative correlations      #
#-------------------------------------------------#
write.csv(df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_7.csv', row.names=FALSE)

write.csv(exp_observed, '~/Desktop/covid_paper_outputs_oasis/Table_S7_exp_observed.csv', row.names=FALSE)
write.csv(highly_predicted, '~/Desktop/covid_paper_outputs_oasis/Table_S7_highlypred.csv', row.names=FALSE)


#----------------------------#
# CROSS CORRELATIONS PLOTS   #
#----------------------------#
setwd('~/Desktop/miRNA/rnaseq/outputs/cross_corr_point_plots')

count = 1

for (i in 1:nrow(df)){
  
  mirna = df[i, 'miRNA']
  mirna = str_replace_all(mirna, "-", '.')
  transcript = df[i, 'Transcript ID']
  gene = df[i, 'Gene']
  
  tmp_df = data[,c(mirna, transcript)]
  names(tmp_df) = c('mirna', 'transcript')
  
  min_x = min(na.omit(tmp_df$mirna))
  min_y = min(na.omit(tmp_df$transcript))
  max_y = max(na.omit(tmp_df$transcript))
  
  p = ggplot(tmp_df, aes(x=mirna, y=transcript)) +
    geom_point(size=2.5) +
    geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_cor(method = "spearman", label.x = min_x, 
             label.y = min_y-0.25,
             size=5.5, color='blue') +
    xlab('miRNA expression') + 
    ylab("Transcript expression") +
    theme_classic() +
    theme(axis.title = element_text(size=16),
          axis.text = element_text(size=16),
          plot.margin = margin(0.5, 1, 0.5, 0.5, 'cm'),
          plot.title = element_text(size=18, hjust=0.5, face='bold')) +
    ggtitle(str_interp("${mirna} - ${gene}")) +
    ylim(c(min_y-0.25, max_y))
  
  assign(str_interp("p${count}"), p)
  
  png(str_interp("${mirna}_${transcript}_${gene}_oasis.png"), width=400, height=300)
  plot(p)
  dev.off()
  
  count = count+1
  
}

p = ggarrange(ggarrange(p1, p3, p6, p9, p13, p26, nrow=2, ncol=3,
                    label.x = 0, label.y = 1, 
                    font.label = list(size=20, color='black', face='bold'),
                    labels = c('A', 'B', 'C', 'D', 'E', 'F')),
              ggarrange(p38, p44, p50, p62, p71, p104 , nrow=2, ncol=3,
                        label.x = 0, label.y = 1, 
                        font.label = list(size=20, color='black', face='bold'),
                        labels = c('G', 'H', 'I', 'J', 'K', 'L')),
              nrow=2, ncol=1, heights = c(0.5, 0.5))


#---------------------------------------------#
# SUPPLEMENTARY FIGURE 10:                    #
# Examples of miRNA-mRNA cross correlations   #
#---------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_10.png', width=900, height=1000)
plot(p)
dev.off()


#------------------------#
# PROCESS IPA RESULTS    #
#------------------------#

# Canonical pathways 
pathways_df = process_canonical_pathways(
  '~/Desktop/miRNA/ipa_results/oasis/blood_pheno/EO_and_HP_canonical_pathways_updated.txt')

write.csv(pathways_df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_8A.csv', 
          row.names=FALSE)


# Diseases & Processes
diseases_df = process_diseases_processes(
  '~/Desktop/miRNA/ipa_results/oasis/blood_pheno/EO_and_HP_diseases_and_functions_updated.txt')

write.csv(diseases_df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_8B.csv', 
          row.names=FALSE)

