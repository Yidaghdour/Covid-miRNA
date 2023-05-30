# Association analyses for all miRNA

rm(list=ls())

#--------------#
# LIBRARIES    #
#--------------#
library(wesanderson)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(dplyr)
library(stringr)
library(MetBrewer)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)

#--------------#
# FUNCTIONS    #
#--------------#

# Functions for binary dependent variables
run_ols_logit = function(y, list_x){
  
  results_df = data.frame(matrix(nrow=length(list_x), ncol=6))
  names(results_df) = c('y', 'x', 'effect', 'p', 'n', 'aic')
  
  k=1
  for (mirna in list_x){
  
    formula = paste(y, "~", mirna)
    res = glm(formula = formula, data=data, family='binomial')

    effect = as.numeric(summary(res)$coefficients[2,1])
    p = as.numeric(summary(res)$coefficients[2,4])
    n = nobs(res)
    aic = as.numeric(summary(res)$aic)
    
    results_df[k, ] = c(y, mirna, effect, p, n, aic)
    k=k+1
    
  }
  
  return(results_df)
}

run_ols_logit_covars = function(y, list_x, covars){
  
  results_df = data.frame(matrix(nrow=length(list_x), ncol=7))
  names(results_df) = c('y', 'x', 'covars', 'effect', 'p', 'n', 'aic')
  
  covars_txt = paste(covars, collapse="+")
  
  k=1
  for (mirna in list_x){
    
    formula = paste(y, "~", mirna)
    
    for (covar in covars){
      formula = paste(formula, "+", covar)
    }
    
    res = glm(formula = formula, data=data, family='binomial')
    
    effect = as.numeric(summary(res)$coefficients[2,1])
    p = as.numeric(summary(res)$coefficients[2,4])
    n = nobs(res)
    aic = as.numeric(summary(res)$aic)
    
    results_df[k, ] = c(y, mirna, covars_txt, effect, p, n, aic)
    k=k+1
    
  }
  
  return(results_df)
}

return_and_save_table_logit = function(results_df, sig, filename){
  
  # Round up numerical columns 
  results_df$p = round(as.numeric(results_df$p), 6)
  results_df$aic = round(as.numeric(results_df$aic), 3)

  # Subset significant results 
  results_df_sig = results_df[results_df$p < sig, ]

  # Order by reducing p-value 
  results_df_sig = results_df_sig[order(results_df_sig$p), ]
  
  # Save resulting table
  write.csv(results_df_sig, str_interp("${filename}.csv"))
  write.csv(results_df_sig, str_interp("${filename}_sig_at_${sig}.csv"))
  
  # Print sumstats
  n_05 = nrow(results_df[results_df$p < 0.05, ])
  n_01 = nrow(results_df[results_df$p < 0.01, ])
  n_001 = nrow(results_df[results_df$p < 0.001, ])
  n_0001 = nrow(results_df[results_df$p < 0.0001, ])
  
  print(str_interp(
    "There are ${n_05} miRNA with p < 0.05, ${n_01} with p < 0.01, ${n_001} with p < 0.001, and ${n_0001} with p < 0.0001."))
    
  return(results_df_sig)
}

save_violins = function(results_df, sig, y_name){
  
  y = unique(results_df$y)
  mirna_to_plot = results_df[results_df$p < sig, 'x']
  
  for (mirna in mirna_to_plot){

    tmp = data[,c(mirna, y)]
    names(tmp) = c('miRNA', 'outcome')
    tmp = tmp[!is.na(tmp$outcome), ]
    
    my_formula = miRNA ~ outcome
    
    p = ggplot(tmp, aes(x=factor(outcome), y=miRNA, fill=factor(outcome))) +
      geom_violin() + 
      stat_summary(fun = "mean",
                   geom = "crossbar", 
                   width = 0.25,
                   colour = "black") +
      ggtitle(mirna) + 
      theme_classic() + ylab('') +
      theme(plot.title = element_text(size=14, vjust=0.5, hjust=0.5, face='bold'),
            legend.position = 'none',
            axis.title = element_text(size=14),
            axis.text = element_text(size=14)) +
      scale_fill_manual(values = c('#E4D1B9', '#A97155')) +
      xlab(y_name)
    
    # Construct filename
    filename = paste(mirna, '_at_', sig, '.png', sep='')
    
    # Save plot
    png(filename, width=500, height=300)
    plot(p)
    dev.off()
    
  }
}

plot_volcano = function(results_df, sig_h1, sig_h2, sig_l, dep_var, save, label=TRUE, hq=FALSE){
  
  results_df = results_df %>%
    
    # Make variables numeric 
    mutate( p = as.numeric(p),
            effect = as.numeric(effect)) %>%
    
    # Annotate miRNA's to highlight 
    mutate( is_annotate = ifelse(p < sig_l, "yes", "no"),
            is_highlight = ifelse(p < sig_h1, "yes", "no"),
            is_highlight_and_label = ifelse(p < sig_h2, 'yes', 'no'))
  
  results_df$log10p = -log10(results_df$p) 
  max_y = max(results_df$log10p)
  
  if (max_y < 3){
    max_y = 3
  }
  
  max_effect = max(abs(results_df$effect))
  
  if (hq == FALSE){
    
    p = ggplot(results_df, aes(x=effect, y=-log10(p))) +
      geom_point(size=2.5) +
      theme_classic() +
      theme(axis.text = element_text(size=13),
            axis.title = element_text(size=13), 
            legend.position = 'none',
            panel.spacing.x = unit(0, "lines"),
            plot.margin = margin(0.5, 0.5, 0.2, 0.5, "cm"),
            plot.title = element_text(size=15, hjust=0.5, face='bold')) +
      xlab('Effect size') + ylab(expression("-"~log[10]~"P")) +
      scale_color_manual(values = colors) +
      ggtitle(dep_var) +
      geom_point(data=subset(results_df, is_highlight=="yes"), 
                 color="#3FA796", size=2.5) +
      geom_point(data=subset(results_df, is_highlight_and_label=="yes"), 
                 color="#C74B50", size=2.5) +
      geom_vline(xintercept = 0, color='grey', linetype='dashed') +
      ylim(c(0, max_y)) +
      xlim(c(-max_effect, max_effect)) + 
      if (label == TRUE){
        geom_label_repel(data=subset(results_df, is_annotate=="yes"), 
                         color="black", aes(label=x), size=4,
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50',
                         max.overlaps = Inf,
                         min.segment.length = 0)
      }
    
  } else {
    p = ggplot(results_df, aes(x=effect, y=-log10(p))) +
      geom_point(size=3.5) +
      theme_classic() +
      theme(axis.text = element_text(size=20),
            axis.title = element_text(size=20), 
            legend.position = 'none',
            panel.spacing.x = unit(0, "lines"),
            plot.margin = margin(0.5, 0.5, 0.2, 0.5, "cm"),
            plot.title = element_text(size=22, hjust=0.5, face='bold')) +
      xlab('Effect size') + ylab(expression("-"~log[10]~"P")) +
      scale_color_manual(values = colors) +
      ggtitle(dep_var) +
      geom_point(data=subset(results_df, is_highlight=="yes"), 
                 color="#3FA796", size=4) +
      geom_point(data=subset(results_df, is_highlight_and_label=="yes"), 
                 color="#C74B50", size=4) +
      geom_vline(xintercept = 0, color='grey', linetype='dashed') +
      ylim(c(0, max_y)) +
      xlim(c(-max_effect, max_effect)) + 
      if (label == TRUE){
        geom_label_repel(data=subset(results_df, is_annotate=="yes"), 
                         color="black", aes(label=x), size=7,
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50',
                         max.overlaps = Inf,
                         min.segment.length = 0)
      }
  }
  
  if (save != 'no'){
    png(save, width=500, height=400)
    plot(p)
    dev.off()
  }
  
  return(p)
}

add_sig_stars = function(table_df, p_col_name){
  
  for (i in 1:nrow(table_df)){
    if (table_df[i, p_col_name] < 0.0001){
      table_df[i, 'significance'] = '****'
    } else if (table_df[i, p_col_name] < 0.001){
      table_df[i, 'significance'] = '***'
    } else if (table_df[i, p_col_name] < 0.01){
      table_df[i, 'significance'] = '**'
    } else if (table_df[i, p_col_name] < 0.05){
      table_df[i, 'significance'] = '*'
    } else if (table_df[i, p_col_name] > 0.05){
      table_df[i, 'significance'] = 'ns'
    }
  }
  
  return(table_df)
  
}

make_replication_plot = function(results_no_covar, results_covar, var, xlab, ylab){
  
  results = combine_sig_results_df(results_no_covar, results_covar)
  
  min = min(min(na.omit(results$effect.x)), min(na.omit(results$effect.y)))
  max = max(max(na.omit(results$effect.x)), max(na.omit(results$effect.y)))
  
  total_max = max(abs(min), abs(max))
  total_tick = floor(total_max)
  
  p = ggplot(results, aes(x=effect.x, y=effect.y, fill=Replicated)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
    geom_hline(yintercept=0, linetype='dashed', color='grey') +
    geom_vline(xintercept=0, linetype='dashed', color='grey') +
    geom_point(shape=21, color='black', size=3) +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          plot.title = element_text(size=17, face='bold', hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=13),
          legend.position = 'top',
          plot.margin = margin(0.5, 0.5, 0.2, 0.5, "cm")) +
    xlab(xlab) +
    ylab(ylab) +
    scale_fill_manual(values = c('#F0F0F0', '#CC704B')) +
    ggtitle("Replication of miRNA associations\nbetween two models") +
    scale_x_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick)) +
    scale_y_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick))
  
  return(p)
}


# Calculate total R2 explained by miRNAs 
calculate_cor_mirnas = function(mirna_list, data, max_cor){
  
  # Subset relevant miRNAs 
  data_subset = data %>%
    select(all_of(mirna_list))
  
  cor_df = as.data.frame(as.table(cor(data_subset, use = 'pairwise.complete.obs')))
  names(cor_df) = c('mirna1', 'mirna2', 'cor')
  
  # Identify highly correlated miRNAs (cor > 0.5)
  tmp_high_cor = cor_df %>%
    filter(cor != 1) %>%
    arrange(cor) %>%
    slice(seq(1, nrow(cor_df), by=2)) %>%
    mutate(cor = abs(as.numeric(cor))) %>%
    filter(cor > max_cor)
  
  # Convert factors to characters
  tmp_high_cor$mirna1 = as.character(tmp_high_cor$mirna1)
  tmp_high_cor$mirna2 = as.character(tmp_high_cor$mirna2)
  
  return(tmp_high_cor)
}

calculate_total_r2_based_on_independent_mirna_logistic = function(data, sig_results_df, max_cor, dep_var){
  
  # All the significant miRNAs 
  sig_mirna_all = sig_results_df$x
  
  # Calculate a correlation matrix of the significant hits 
  tmp_high_cor = calculate_cor_mirnas(sig_mirna_all, data, max_cor)
  
  # Make a list for independent miRNAs 
  sig_mirna = sig_mirna_all
  
  # Remove miRNAs to make a list of independent ones
  while (nrow(tmp_high_cor) > 1){
    
    # Make a table of most frequently occuring miRNAs and remove those first 
    tmp_freq = data.frame(table(c(as.character(tmp_high_cor$mirna1), 
                                  as.character(tmp_high_cor$mirna2))))
    tmp_freq = tmp_freq[order(tmp_freq$Freq, decreasing = TRUE),]
    row.names(tmp_freq) = seq(1, nrow(tmp_freq))
    
    # Find the most commonly correlated miRNA 
    mirna = as.character(tmp_freq[1, 'Var1'])
    
    # Remove it from significant miRNAs and redo the whole thing 
    sig_mirna = sig_mirna[-which(sig_mirna == mirna)]
    
    # Recalculate correlations and all 
    tmp_high_cor = calculate_cor_mirnas(sig_mirna, data, max_cor)
    
  }
  
  # Calculate R2 based on independent miRNAs 
  mirnas_txt = paste(sig_mirna, collapse=" + ")
  full_formula = str_interp("${dep_var} ~ Age_numeric + Gender_binary + onset_to_admit + ${mirnas_txt}")
  
  results = glm(full_formula, data, na.action = 'na.omit', family='binomial')
  results_age_sex = glm(str_interp("${dep_var} ~ Age_numeric + Gender_binary"), na.action = 'na.omit', data, family='binomial')
  results_null_model = glm(str_interp("${dep_var} ~ 1"), na.action = 'na.omit', data, family='binomial')
  
  macfadden_r2 = round(1-(as.numeric(logLik(results)))/as.numeric(logLik(results_null_model)), 4)
  results_n = length(sig_mirna)
  
  macfadden_r2_over_gender_sex = round(1-(as.numeric(logLik(results)))/as.numeric(logLik(results_age_sex)), 4)
  
  print(str_interp("The ${results_n} independent miRNAs (cor < ${max_cor}) + Age, Sex and Onset-to-admit explain ${macfadden_r2*100}% of the variation in ${dep_var}."))
  print(str_interp("The ${results_n} independent miRNAs (cor < ${max_cor}) + explain ${macfadden_r2_over_gender_sex*100}% of the variation in ${dep_var} over age, sex and onset to admit only."))
  print(str_interp("The miRNAs are ${sig_mirna}."))  
  
  return(results)
}

combine_sig_results_df = function(df_nocovar, df_covar){
  
  df_nocovar$covars = 'None'
  
  merged_results = merge(df_nocovar, df_covar, by='x', all=TRUE)
  merged_results[is.na(merged_results)] = 0
  
  merged_results$Replicated = 'Replicated'
  merged_results[merged_results$effect.x == 0 | merged_results$effect.y == 0, 'Replicated'] = 'Not replicated'
  
  merged_results$effect.x = as.numeric(merged_results$effect.x)
  merged_results$effect.y = as.numeric(merged_results$effect.y)
  
  return(merged_results)
  
}

add_fdr_correction = function(df, p_colname, correction_method){
  
  df$FDR_P = p.adjust(df[,p_colname], 
                      method=correction_method, n=nrow(df))
  
  n_sig_nominal = nrow(df %>% filter(p < 0.05))
  n_sig_fdr = nrow(df %>% filter(FDR_P < 0.05))
  
  print(str_interp("There are ${n_sig_nominal} associations with nominal P < 0.05."))
  print(str_interp("There are ${n_sig_fdr} associations with FDR P < 0.05."))
  
  return(df)
  
}


#--------#
# LISTS  #
#--------#
colors = met.brewer("Tiepolo", 9)


#--------------#
# LOAD DATA    #
#--------------#
setwd('~/Desktop/miRNA/mirna_processed/n96/')
data = read.csv('mirna_oasis_5_in_50_percent_standard_clinical.csv')

var_names = data.frame(names(data))
mirna_list = names(data)[3:634]


#-------------------------------------------------#
# ICU ADMISSION ANALYSES - NOT gender-corrected   #
#-------------------------------------------------#

# 1. MIRNA vs out_ICU (not gender-corrected) -----------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_icu')

covars = c('Age_numeric', 'onset_to_admit')

# Compute OLS logit results 
results_df = run_ols_logit('out_ICU', mirna_list)
results_covars_df = run_ols_logit_covars('out_ICU', mirna_list, covars)

# Save results as a table 
sig_results_df_icu = return_and_save_table_logit(results_df, 0.05, 'table_ICU')
sig_results_covars_df_icu = return_and_save_table_logit(results_covars_df, 0.05, 'table_ICU')

# Fix up table for plotting
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")
results_covars_df$x = str_replace_all(results_covars_df$x, "hsa-", "")

# Plot results (p<0.05 highlighted, p<0.005 labeled)
p_icu = plot_volcano(results_df, 0.05, 0.01, 0.01, 'ICU admission', 'no')
p_icu_covar = plot_volcano(results_covars_df, 0.05, 0.01, 0.00315, 'miRNAs associated with ICU admission', 'no')
p_icu_covar_hq = plot_volcano(results_covars_df, 0.05, 0.01, 0.00315, 
                              'miRNAs associated with ICU admission', 'no', hq=TRUE)

#---------------------------#
# FIGURE 1C: ICU ADMISSION  #
#---------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Figure_1C.png', width=400, height=380)
plot(p_icu_covar)
dev.off()

png('~/Desktop/covid_paper_outputs_oasis/Figure_1C_hq.png', width=600, height=570)
plot(p_icu_covar_hq)
dev.off()




#-------------------------------------------#
# SUPPLEMENTARY TABLE 4A:                   #
# micro RNAs associated with ICU admission  #
#-------------------------------------------#
sig_results_covars_df_icu_for_table = add_sig_stars(sig_results_covars_df_icu, 'p')
names(sig_results_covars_df_icu_for_table) = c('Dependent variable', 'miRNA', 'Covars', 'Effect', 
                                     'P-value', 'N', 'AIC', 'Significance')

sig_results_covars_df_icu_for_table$miRNA = str_replace_all(sig_results_covars_df_icu_for_table$miRNA,
                                                  '\\.', "-")

write.csv(sig_results_covars_df_icu_for_table, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4A.csv')



#---------------------------------------------#
# ICU ADMISSION ANALYSES - gender-corrected   #
#---------------------------------------------#
setwd('~/Desktop/mirna/gwas_oasis/mirna_icu')

covars_w_gender = c('Age_numeric', 'Gender_binary', 'onset_to_admit')

# Compute OLS logit results 
results_covars_w_gender_df = run_ols_logit_covars('out_ICU', mirna_list, covars_w_gender)

# Save results as a table 
sig_results_covars_w_gender_df_icu = 
  return_and_save_table_logit(results_covars_w_gender_df, 0.05, 'table_ICU')

# Fix up table for plotting
results_covars_w_gender_df$x = str_replace_all(results_covars_w_gender_df$x, "\\.", "-")
results_covars_w_gender_df$x = str_replace_all(results_covars_w_gender_df$x, "hsa-", "")

# Plot results (p<0.05 highlighted, p<0.005 labeled)
p_icu_w_gender_covar = plot_volcano(results_covars_w_gender_df, 
                                     0.05, 0.01, 0.003, 
                                     'miRNAs associated with ICU admission\n(adjusted for gender)', 'no')


## Compare with and without gender as a covariate 
p_rep_icu = make_replication_plot(sig_results_covars_df_icu, 
                                  sig_results_covars_w_gender_df_icu, 'ICU',
                                  'Effect size of miRNA, not corrected for gender',
                                  'Effect size of miRNA, corrected for gender')

p = ggarrange(p_icu_w_gender_covar, p_rep_icu, nrow=1,
              label.x = 0, label.y = 1, 
              font.label = list(size=20, color='black', face='bold'),
              labels = c('A', 'B'))


#--------------------------------------------------------#
# SUPPLEMENTARY FIGURE 3A-B:                             #
# (A) miRNA associations in a model adjusted for gender  #
# (B) Replication with and without gender as covariate   #
#--------------------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_3.png', width=700, height=370)
plot(p)
dev.off()


#-------------------------------------------#
# SUPPLEMENTARY TABLE 4B:                   #
# micro RNAs associated with ICU admission  #
#-------------------------------------------#
sig_results_covars_w_gender_df_icu_for_table = add_sig_stars(sig_results_covars_w_gender_df_icu, 'p')
names(sig_results_covars_w_gender_df_icu_for_table) = c('Dependent variable', 'miRNA', 'Covars', 'Effect', 
                                     'P-value', 'N', 'AIC', 'Significance')

sig_results_covars_w_gender_df_icu_for_table$miRNA = str_replace_all(sig_results_covars_w_gender_df_icu_for_table$miRNA,
                                                  '\\.', "-")

write.csv(sig_results_covars_w_gender_df_icu_for_table, 
          '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4B.csv')


#-----------------------------------------------------------#
# IN-TEXT:                                                  #
# Overlap of miRNAs between gender- and non-gender model    #
#-----------------------------------------------------------#

sig_p01_not_gender_corrected = sig_results_covars_df_icu %>%
  filter(p < 0.01)

sig_p01_w_gender_corrected = sig_results_covars_w_gender_df_icu %>%
  filter(p < 0.01)

shared_mirna =
  intersect(sig_p01_not_gender_corrected$x, sig_p01_w_gender_corrected$x)

sig_p01_not_gender_corrected$is_shared = FALSE
sig_p01_not_gender_corrected[which(sig_p01_not_gender_corrected$x %in% shared_mirna),
                            'is_shared'] = TRUE 

sig_p01_w_gender_corrected$is_shared = FALSE
sig_p01_w_gender_corrected[which(sig_p01_not_gender_corrected$x %in% shared_mirna),
                          'is_shared'] = TRUE 


