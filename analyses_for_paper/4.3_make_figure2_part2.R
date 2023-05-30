rm(list=ls())

library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(dplyr)
library(stringr)
library(ggpubr)


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

run_ols_covars = function(y, list_x, rm_outliers=TRUE, n_sd=3, covars, top_mirna){
  
  results_df = data.frame(matrix(nrow=length(list_x), ncol=7))
  names(results_df) = c('y', 'x', 'covars', 'effect', 'p', 'n', 'r2_adjusted')
  
  covars_txt = paste(covars, collapse="+")
  
  if (rm_outliers == TRUE){
    
    mean_y = mean(na.omit(data[,y]))
    sd_y = sd(na.omit(data[,y]))
    
    # Remove data points that are 2 SD away from mean (outside of middle 95%)
    tmp_data = data[data[,y] < mean_y+n_sd*sd_y,]
    tmp_data = tmp_data[tmp_data[,y] > mean_y-n_sd*sd_y,]
    
    n_removed = nrow(data) - nrow(tmp_data)
    print(str_interp("${n_removed} datapoints were removed as outliers."))
  }
  
  k=1
  for (mirna in list_x){
    
    formula = paste(y, "~", mirna)
    
    for (covar in covars){
      formula = paste(formula, "+", covar)
    }
    
    if (rm_outliers == TRUE){
      res = lm(formula = formula, data=tmp_data)
    } else {
      res = lm(formula = formula, data=data)
    }
    
    beta = as.numeric(summary(res)$coefficients[2,1])
    p = as.numeric(summary(res)$coefficients[2,4])
    n = nobs(res)
    r2_adj = as.numeric(summary(res)$adj.r.squared)
    
    results_df[k, ] = c(y, mirna, covars_txt, beta, p, n, r2_adj)
    k=k+1
    
  }
  
  # Find top miRNA (lowest P)
  results_df$p = as.numeric(results_df$p)

  # Annotate top miRNA 
  results_df$top_mirna = 'no'
  
  if (min(results_df$p) < 0.01){
    results_df[which(results_df$x == top_mirna), 'top_mirna'] = 'yes'
  }
  
  return(results_df)
}

plot_volcano = function(results_df, sig_h1, sig_h2, dep_var, save, label=TRUE){
  
  results_df = results_df %>%
    
    # Make variables numeric 
    mutate( p = as.numeric(p),
            effect = as.numeric(effect)) %>%
    
    # Annotate miRNA's to highlight 
    mutate( is_annotate = ifelse(p < sig_h2, "yes", "no"),
            is_highlight = ifelse(p < sig_h1, "yes", "no"),
            is_labeled = ifelse(top_mirna == 'yes', 'yes', 'no'))
  
  results_df$log10p = -log10(results_df$p) 
  max_y = max(results_df$log10p)
  
  if (max_y < 3){
    max_y = 3
  }
  
  max_effect = max(abs(results_df$effect))
  
  p = ggplot(results_df, aes(x=effect, y=-log10(p))) +
    geom_point(size=2.5) +
    theme_classic() +
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=16), 
          legend.position = 'none',
          plot.margin = margin(0.5, 1, 0.5, 0.5, 'cm'),
          plot.title = element_text(size=18, hjust=0.5, face='bold')) +
    xlab('Effect size') + ylab(expression("-"~log[10]~"P")) +
    scale_color_manual(values = colors) +
    ggtitle(str_interp("${dep_var}")) +
    geom_point(data=subset(results_df, is_highlight=="yes"), 
               color="#3FA796", size=2.5) +
    geom_point(data=subset(results_df, is_annotate=="yes"), 
               color="#C74B50", size=2.5) +
    geom_vline(xintercept = 0, color='grey', linetype='dashed') +
    ylim(c(0, max_y)) +
    xlim(c(-max_effect, max_effect)) + 
    if (label == TRUE){
      geom_label_repel(data=subset(results_df, is_labeled=="yes"), 
                       color="black", aes(label=x), size=5.5,
                       box.padding   = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'grey50',
                       max.overlaps = Inf,
                       min.segment.length = 0)
    }
  
  
  
  if (save != 'no'){
    png(save, width=500, height=400)
    plot(p)
    dev.off()
  }
  
  return(p)
}

return_and_save_table = function(results_df, sig, filename){
  
  # Round up numerical columns 
  results_df$effect = round(as.numeric(results_df$effect), 4)
  results_df$p = round(as.numeric(results_df$p), 6)
  results_df$r2_adjusted = round(as.numeric(results_df$r2_adjusted), 3)
  
  # Subset significant results 
  results_df_sig = results_df[results_df$p < sig, ]
  
  # Order by reducing p-value 
  results_df_sig = results_df_sig[order(results_df_sig$p), ]
  
  # Print sumstats
  n_05 = nrow(results_df[results_df$p < 0.05, ])
  n_01 = nrow(results_df[results_df$p < 0.01, ])
  n_001 = nrow(results_df[results_df$p < 0.001, ])
  n_0001 = nrow(results_df[results_df$p < 0.0001, ])
  
  print(str_interp(
    "There are ${n_05} miRNA with p < 0.05, ${n_01} with p < 0.01, ${n_001} with p < 0.001, and ${n_0001} with p < 0.0001."))
  
  # Save resulting table
  write.csv(results_df_sig, str_interp("${filename}.csv"))
  write.csv(results_df_sig, str_interp("${filename}_sig_at_${sig}.csv"))
  
  return(results_df_sig)
}


#-------------#
# LOAD DATA   #
#-------------#
data = read.csv('~/Desktop/miRNA/rnaseq/all_data_mirna_rna_oasis.csv')
data = standardize_blood_phenotypes(data)

rnaseq_ann_raw = read.delim2('~/Desktop/miRNA/rnaseq/TPM/NGS007_T1_S198.tpm')
rnaseq_ann = rnaseq_ann_raw %>%
  select(transcript_id, ref_gene_name) %>%
  mutate(ref_gene_name = trimws(ref_gene_name),
         transcript_id = trimws(transcript_id)) 

med_results = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_7.csv')


#------------#
# VARIABLES  #
#------------#
mirna_list = names(data)[3:634]
covars = c('Age_numeric', 'Gender_binary')

bcl2_transcripts = (rnaseq_ann %>%
                      filter(ref_gene_name == 'BCL2'))$transcript_id

top_mirna = 'hsa.miR.143.3p'


#------------#
# ANALYSIS   #
#------------#

# Run OLS 
results_covars_df = run_ols_covars('neutrophil_scaled', mirna_list, TRUE, 3, covars, top_mirna)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_covars_df_neutrophil = return_and_save_table(results_covars_df, 0.05, 
                                                         'table_neutrophil_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_neutrophil_covars = plot_volcano(results_covars_df, 0.05, 0.01, 
                                   'Absolute neutrophil count', 'no', TRUE)


# Ploting correlations
mirna_txt = med_results[1,'miRNA']
mirna = str_replace_all(mirna_txt, "-", ".")

transcript = med_results[1,'Transcript.ID']
gene = med_results[1,'Gene']

pheno = 'neutrophil_scaled'
pheno_txt = 'Neutrophil count'

tmp_df = data[,c(mirna, transcript, pheno)]
names(tmp_df) = c('mirna', 'transcript', 'pheno')

min_x = min(na.omit(tmp_df$mirna))
max_x = max(na.omit(tmp_df$mirna))

min_y = min(na.omit(tmp_df$transcript))
max_y = max(na.omit(tmp_df$transcript))

p_mirna_gene = ggplot(tmp_df, aes(x=mirna, y=transcript)) +
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
  ggtitle(str_interp("${mirna_txt} - ${gene}")) +
  ylim(c(min_y-0.25, max_y))

p_gene_pheno = ggplot(tmp_df, aes(x=transcript, y=pheno)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm', formula= y~x, se=FALSE) +
  stat_cor(method = "spearman", label.x = min_y, 
           label.y = min_y-1,
           size=5.5, color='blue') +
  xlab('Transcript expression') + 
  ylab("Neutrophil count") +
  theme_classic() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, 'cm'),
        plot.title = element_text(size=18, hjust=0.5, face='bold')) +
  ggtitle(str_interp("${gene} - ${pheno_txt}")) +
  ylim(c(min_y-1, max_y))


#----------------------------------------------#
# FIGURE 2                                     #
# miRNA regulating neutrophil count via BCL2   #
#----------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Figure_2.png', width=1050, height=300)
ggarrange(p_neutrophil_covars, p_mirna_gene, p_gene_pheno,
          ncol=3, nrow = 1, widths = c(0.3, 0.35, 0.35),
          label.x = 0, label.y = 1, labels = c('A', 'B', 'C'),
          font.label = list(size=22, color='black', face='bold'))
dev.off()


