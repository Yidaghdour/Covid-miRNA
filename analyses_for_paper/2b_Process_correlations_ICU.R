rm(list=ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)

#---------#
# LISTS   #
#---------#
ipa_results_dir = '~/Desktop/miRNA/ipa_results/oasis/'


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
load_icu_ipa_data = function(ipa_results_dir, filename){
  
  setwd(ipa_results_dir)
  
  raw_main_df = readxl::read_xls(str_interp(filename))
  raw_main_df = raw_main_df[-1,]
  
  names(raw_main_df) = c('mirna_id', 'mirna_id_explained', 
                         'source', 'confidence', 'gene', 'pathway')
  
  main_df = raw_main_df %>%
    select(mirna_id, confidence, gene, pathway) %>%
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


#-------------#
# LOAD DATA   #
#-------------#

# Load the correlation results
icu_s_corrs = read.csv('~/Desktop/miRNA/rnaseq/outputs/icu_corrs_spearman_sig_fdr05_oasis.csv')
icu_s_corrs$mirna = str_replace_all(icu_s_corrs$mirna, '\\.', '-')

icu_mirna = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4A.csv')
icu_mirna = icu_mirna %>%
  filter(P.value < 0.01)

data = read.csv('~/Desktop/miRNA/rnaseq/all_data_mirna_rna_oasis.csv')


#--------------------------#
# ANNOTATE miRNA TARGETS   #
#--------------------------#

# IPA targets
main_df = load_icu_ipa_data(ipa_results_dir, 'icu_mirna_targets_updated.xls')

# Annotate miRNA-RNA pairs 
icu_s_corrs_ann = annotate_mirna_table(icu_s_corrs)

write.csv(icu_s_corrs_ann, 
          '~/Desktop/miRNA/rnaseq/outputs/corrs_spearman_sig_fdr05_icu_oasis_annotated_updated.csv')


#---------------------------------#
# PROCESS DATA (informs ms text)  #
#---------------------------------#

## NEGATIVELY CORRELATED miRNA-mRNA 
icu_s_corrs_neg_ann = icu_s_corrs_ann %>%
  filter(corr < 0)

icu_s_corrs_neg_ann = add_sig_stars(icu_s_corrs_neg_ann, 'FDR_P', 'Significance')

n_corrs = nrow(icu_s_corrs_ann)
n_negative_corrs = nrow(icu_s_corrs_neg_ann)
min_neg_cor = round(min(icu_s_corrs_neg_ann$corr), 2)
max_neg_cor = round(max(icu_s_corrs_neg_ann$corr), 2)
mean_neg_cor = round(mean(icu_s_corrs_neg_ann$corr), 2)
sd_neg_cor = round(sd(icu_s_corrs_neg_ann$corr), 2)

print(str_interp(
  "There are a total of ${n_corrs} correlated miRNA-mRNA pairs."))

print(str_interp(
  "There are a total of ${n_negative_corrs} negatively-correlated miRNA-mRNA pairs."))

print(str_interp(
  "The range was ${min_neg_cor} to ${max_neg_cor}. The mean was ${mean_neg_cor} (SD = ${sd_neg_cor})"))

n_unique_mirna = length(unique(icu_s_corrs_neg_ann$mirna))

print(str_interp("It concerned ${n_unique_mirna} miRNAs."))


## HOW MANY miRNAS HAVE TARGETS FROM IPA? 
n_unique_mirna = length(unique(icu_s_corrs_neg_ann$mirna))

tmp = icu_s_corrs_neg_ann %>%
  filter(predicted_target != 'No information')

n_mirna_w_target_data = length(unique(tmp$mirna))

print(str_interp(
  "We have target data for ${n_mirna_w_target_data} out of ${n_unique_mirna} miRNAs."))


## EXPERIMENTALLY OBSERVED 
exp_observed = icu_s_corrs_neg_ann %>%
  filter(confidence == 'Experimentally Observed')

n_exp_obs = nrow(exp_observed)
n_mirna = length(unique(exp_observed$mirna))
n_genes = length(unique(exp_observed$ref_gene_name))

print(str_interp(
  "There were ${n_exp_obs} experimentally observed pairs, concerning ${n_mirna} miRNAs and ${n_genes} genes."))


## HIGHLY PREDICTED 
highly_predicted = icu_s_corrs_neg_ann %>%
  filter(confidence == 'High (predicted)') %>%
  select(mirna, ref_gene_name) %>%
  distinct()

n_exp_obs = nrow(highly_predicted)
n_mirna = length(unique(highly_predicted$mirna))
n_genes = length(unique(highly_predicted$ref_gene_name))

print(str_interp(
  "There were ${n_exp_obs} highly predicted pairs, concerning ${n_mirna} miRNAs and ${n_genes} genes."))


## EXPERIMENTALLY OBSERVED + HIGHLY PREDICTED 
df = icu_s_corrs_neg_ann %>%
  filter(confidence %in% c('Experimentally Observed',
                           'High (predicted)')) %>%
  select(mirna, rna, ref_gene_name, corr, FDR_P, Significance, confidence) %>%
  arrange(confidence, mirna)

names(df) = c('miRNA', 'Transcript ID', 'Gene', 'Correlation coefficient',
              'FDR P', 'Significance', 'Confidence')

df$`Correlation coefficient` = round(df$`Correlation coefficient`, 3)


#-------------------------------------------------#
# SUPPLEMENTARY TABLE 5A:                         #
# Experimentally observed and highly predicted    #
# miRNA-RNA pairs with negative correlations      #
#-------------------------------------------------#
write.csv(df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_5A.csv', row.names=FALSE)


#----------------------------#
# CROSS CORRELATIONS PLOTS   #
#----------------------------#
setwd('~/Desktop/miRNA/rnaseq/outputs/cross_corr_point_plots')

count = 1

for (i in 1:nrow(df)){
  
  mirna = df[i, 'miRNA']
  mirna = str_replace_all(mirna, '-', '.')
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

p = ggarrange(p1, p2, p9, p10, p11, p12, nrow=2, ncol=3,
              label.x = 0, label.y = 1, 
              font.label = list(size=20, color='black', face='bold'),
              labels = c('A', 'B', 'C', 'D', 'E'))


#---------------------------------------------#
# SUPPLEMENTARY FIGURE 7:                     #
# Examples of miRNA-mRNA cross correlations   #
#---------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_7.png', width=900, height=500)
plot(p)
dev.off()
