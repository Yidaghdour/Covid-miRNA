# Association analyses for all miRNA

rm(list=ls())

#--------------#
# LIBRARIES    #
#--------------#
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(dplyr)
library(stringr)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)


#--------------#
# FUNCTIONS    #
#--------------#

# Functions for continuous dependent variables
run_ols = function(y, list_x, rm_outliers=TRUE, n_sd=3){
  
  results_df = data.frame(matrix(nrow=length(list_x), ncol=6))
  names(results_df) = c('y', 'x', 'effect', 'p', 'n', 'r2_adjusted')
  
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
    
    if (rm_outliers == TRUE){
      res = lm(formula = formula, data=tmp_data)
    } else {
      res = lm(formula = formula, data=data)
    }
    
    beta = as.numeric(summary(res)$coefficients[2,1])
    p = as.numeric(summary(res)$coefficients[2,4])
    n = nobs(res)
    r2_adj = as.numeric(summary(res)$adj.r.squared)
    
    results_df[k, ] = c(y, mirna, beta, p, n, r2_adj)
    k=k+1
    
  }
  
  # Find top miRNA (lowest P)
  results_df$p = as.numeric(results_df$p)
  top_mirna = results_df[results_df$p == min(results_df$p), 'x']
  
  # Annotate top miRNA 
  results_df$top_mirna = 'no'
  
  if (min(results_df$p) < 0.01){
    results_df[which(results_df$x == top_mirna), 'top_mirna'] = 'yes'
  }
  
  return(results_df)
}

run_ols_covars = function(y, list_x, rm_outliers=TRUE, n_sd=3, covars){
  
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
  top_mirna = results_df[results_df$p == min(results_df$p), 'x']
  
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
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13), 
          legend.position = 'none',
          panel.spacing.x = unit(0, "lines"),
          plot.title = element_text(size=15, hjust=0.5, face='bold')) +
    xlab('Effect size') + ylab(expression("-"~log[10]~"P")) +
    scale_color_manual(values = colors) +
    ggtitle(str_interp("\n${dep_var}")) +
    geom_point(data=subset(results_df, is_highlight=="yes"), 
               color="#3FA796", size=2.5) +
    geom_point(data=subset(results_df, is_annotate=="yes"), 
               color="#C74B50", size=2.5) +
    geom_vline(xintercept = 0, color='grey', linetype='dashed') +
    ylim(c(0, max_y)) +
    xlim(c(-max_effect, max_effect)) + 
    if (label == TRUE){
      geom_label_repel(data=subset(results_df, is_labeled=="yes"), 
                       color="black", aes(label=x), size=4,
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

save_scatters = function(results_df, sig, y_name, rm_outliers, n_sd){
  
  y = unique(results_df$y)
  mirna_to_plot = results_df[as.numeric(results_df$p) < sig, 'x']
  
  for (mirna in mirna_to_plot){
    
    if (rm_outliers == TRUE){
      
      mean_y = mean(na.omit(data[,y]))
      sd_y = sd(na.omit(data[,y]))
      
      # Remove data points that are n_sd SD away from mean (outside of middle 95%)
      tmp_data = data[data[,y] < mean_y+n_sd*sd_y,]
      tmp_data = tmp_data[tmp_data[,y] > mean_y-n_sd*sd_y,]
      
      n_removed = nrow(data) - nrow(tmp_data)
      #print(str_interp("${n_removed} datapoints were removed as outliers for ${y_name}."))
      
      tmp = tmp_data[,c(y, mirna)]
      
    } else {
      tmp = data[,c(y, mirna)]
    }
    
    names(tmp) = c('outcome', 'miRNA')
    
    my_formula = outcome ~ miRNA
    
    p = ggplot(tmp, aes(x=miRNA, y=outcome)) +
      geom_point() +
      geom_smooth(formula = y~x, method = "lm", se=TRUE) +
      stat_poly_eq(formula = y~x,
                   aes(label = ..rr.label..), 
                   parse = TRUE, size = 4.5) +
      theme_classic() +
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14),
            plot.title = element_text(size=15, face='bold', hjust=0.5)) +
      ggtitle(mirna) +
      xlab('') + ylab(y_name)
    
    # Construct filename
    filename = paste(mirna, '_at_', sig, '.png', sep='')
    
    # Save plot
    png(filename, width=500, height=300)
    plot(p)
    dev.off()
    
  }
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

calculate_total_r2_based_on_independent_mirna = function(data, sig_results_df, max_cor, dep_var){
  
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
  full_formula = str_interp("${dep_var} ~ Age_numeric + Gender_binary + ${mirnas_txt}")
  age_sex_formula = str_interp("${dep_var} ~ Age_numeric + Gender_binary")
  
  results = summary(lm(full_formula, data, na.action = 'na.omit'))
  results_adj_r2 = round(results$adj.r.squared, 4)
  results_n = length(sig_mirna)
  
  results_age_sex = summary(lm(age_sex_formula, data, na.action = 'na.omit'))
  
  incremental_r2 = round(results$adj.r.squared - results_age_sex$adj.r.squared, 4)
  
  print(str_interp("The ${results_n} independent miRNAs (cor < ${max_cor}) + Age and sex explain ${results_adj_r2*100}% of the variation in ${dep_var}."))
  print(str_interp("The ${results_n} independent miRNAs (cor < ${max_cor}) explain ${incremental_r2*100}% over age and sex."))
  print(str_interp("The miRNAs are ${sig_mirna}."))
  
  return(results)
}


# Replication plots
make_replication_plot = function(results_no_covar, results_covar, var){
  
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
    xlab('Effect size of miRNA, no covariates') +
    ylab('Effect size of miRNA, with covariates') +
    scale_fill_manual(values = c('#F0F0F0', '#CC704B')) +
    ggtitle(str_interp("${var}")) +
    scale_x_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick)) +
    scale_y_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick))
  
  return(p)
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
  
  df$p = as.numeric(df$p)
  
  df$FDR_P = p.adjust(df[,p_colname], 
                      method=correction_method, n=nrow(df))
  
  n_sig_nominal = nrow(df %>% filter(p < 0.05))
  n_sig_fdr = nrow(df %>% filter(FDR_P < 0.1))
  
  print(str_interp("There are ${n_sig_nominal} associations with nominal P < 0.05."))
  print(str_interp("There are ${n_sig_fdr} associations with FDR P < 0.1."))
  
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


#--------------------------------------------------------------------#
# BLOOD BIOMARKERS: Urea, Chloride, CRP, IL6, Neutrophil count,      #
#                   Lymphocyte count, Neutro/Lympho ratio, D-dimers  #
#--------------------------------------------------------------------#
covars = c('Age_numeric', 'Gender_binary')


## UREA -----------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_urea')

# Prepare dependent variable
data$FIRST_UREA_LVL[data$FIRST_UREA_LVL == 0] = NA
data$urea_scaled = scale(data$FIRST_UREA_LVL, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('urea_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('urea_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_urea = return_and_save_table(results_df, 0.05, 'table_urea_scaled')
sig_results_covars_df_urea = return_and_save_table(results_covars_df, 0.05, 
                                                   'table_urea_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_urea = plot_volcano(results_df, 0.05, 0.01, 'Urea', 'no', FALSE)
p_urea_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Urea', 'no', TRUE)


## CHLORIDE -----------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_chloride')

# Prepare dependent variable
data$FIRST_CHLORIDE_LVL[data$FIRST_CHLORIDE_LVL == 0] = NA
data$chloride_scaled = scale(data$FIRST_CHLORIDE_LVL, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('chloride_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('chloride_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_chloride = return_and_save_table(results_df, 0.05, 'table_chloride_scaled')
sig_results_covars_df_chloride = return_and_save_table(results_covars_df, 0.05, 
                                                'table_chloride_covars_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_chloride = plot_volcano(results_df, 0.05, 0.01, 'Chloride', 'no', FALSE)
p_chloride_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Chloride', 
                                 'no', TRUE)


## CRP ------------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_crp')

# Prepare dependent variable
data$FIRST_C_REACTIVE_PROT[data$FIRST_C_REACTIVE_PROT == 0] = NA
data$crp_scaled = scale(data$FIRST_C_REACTIVE_PROT, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('crp_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('crp_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_crp = return_and_save_table(results_df, 0.05, 'table_crp_scaled')
sig_results_covars_df_crp = return_and_save_table(results_covars_df, 0.05, 
                                                  'table_crp_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_crp = plot_volcano(results_df, 0.05, 0.01, 'C Reactive Protein', 'no', FALSE)
p_crp_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'C Reactive Protein', 
                            'no', TRUE)


## IL6 ------------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_il6')

# Prepare dependent variable
data$FIRST_INTERLEUKIN6[data$FIRST_INTERLEUKIN6 == 0] = NA
data$il6_scaled = scale(data$FIRST_INTERLEUKIN6, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('il6_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('il6_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_il6 = return_and_save_table(results_df, 0.05, 'table_il6_scaled')
sig_results_covars_df_il6 = return_and_save_table(results_covars_df, 0.05, 
                                                  'table_il6_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_il6 = plot_volcano(results_df, 0.05, 0.01, 'Interleukin 6', 'no', FALSE)
p_il6_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Interleukin 6', 
                            'no', TRUE)


## NEUTROPHIL COUNT -----------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_neutro')

# Prepare dependent variable
data$FIRST_ABSOLUTE_NEUTROPHIL[data$FIRST_ABSOLUTE_NEUTROPHIL == 0] = NA
data$neutrophil_scaled = scale(data$FIRST_ABSOLUTE_NEUTROPHIL, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('neutrophil_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('neutrophil_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_neutrophil = return_and_save_table(results_df, 0.05, 'table_neutrophil_scaled')
sig_results_covars_df_neutrophil = return_and_save_table(results_covars_df, 0.05, 
                                                         'table_neutrophil_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_neutrophil = plot_volcano(results_df, 0.05, 0.01, 'Absolute neutrophil count', 'no',
                            FALSE)
p_neutrophil_covars = plot_volcano(results_covars_df, 0.05, 0.01, 
                                   'Absolute neutrophil count', 'no', TRUE)


## LYMPHOCYTE COUNT -----------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_lympho')

# Prepare dependent variable
data$FIRST_ABSOLUTE_LYMPHOCYTE[data$FIRST_ABSOLUTE_LYMPHOCYTE == 0] = NA
data$lymphocyte_scaled = scale(data$FIRST_ABSOLUTE_LYMPHOCYTE, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('lymphocyte_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('lymphocyte_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_lymphocyte = return_and_save_table(results_df, 0.05, 'table_lymphocyte_scaled')
sig_results_covars_df_lymphocyte = return_and_save_table(results_covars_df, 0.05, '
                                                         table_lymphocyte_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_lymphocyte = plot_volcano(results_df, 0.05, 0.01, 'Absolute lymphocyte count', 'no',
                            FALSE)
p_lymphocyte_covars = plot_volcano(results_covars_df, 0.05, 0.01, 
                                   'Absolute lymphocyte count', 'no', TRUE)


## NEUTROPHIL-TO-LYMPHOCYTE RATIO -------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_nl_ratio')

# Prepare dependent variable
data$FIRST_NEUTRO_TO_LYMPHO_RATIO[data$FIRST_NEUTRO_TO_LYMPHO_RATIO == 0] = NA
data$neutro_lympho_ratio_scaled = scale(data$FIRST_NEUTRO_TO_LYMPHO_RATIO, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('neutro_lympho_ratio_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('neutro_lympho_ratio_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_nl_ratio = return_and_save_table(results_df, 0.05, 'table_nl_ratio_scaled')
sig_results_covars_df_nl_ratio = return_and_save_table(results_covars_df, 0.05, '
                                                         table_nl_ratio_scaled')

# Plot results (p<0.05 highlighted, p<0.005 labeled)
p_nl_ratio = plot_volcano(results_df, 0.05, 0.01, 'Neutrophil-to-lymphocyte ratio', 'no',
                            FALSE)
p_nl_ratio_covars = plot_volcano(results_covars_df, 0.05, 0.01, 
                                   'Neutrophil-to-lymphocyte ratio', 'no', TRUE)


## D-DIMERS COUNT -------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_ddimer')

# Prepare dependent variable
data$FIRST_D_DIMER[data$FIRST_D_DIMER == 0] = NA
data$ddimer_scaled = scale(data$FIRST_D_DIMER, center=TRUE, scale=TRUE)

# Compute OLS results 
results_df = run_ols('ddimer_scaled', mirna_list, TRUE, 3)
results_covars_df = run_ols_covars('ddimer_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_df_ddimer = return_and_save_table(results_df, 0.05, 'table_ddimer_scaled')
sig_results_covars_df_ddimer = return_and_save_table(results_covars_df, 0.05, 
                                                     'table_ddimer_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_ddimer = plot_volcano(results_df, 0.05, 0.01, 'D-dimers', 'no', FALSE)
p_ddimer_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'D-dimers', 'no', TRUE)


#------------------------------------------------------------------#
# SUPPLEMENTARY FIGURE 6:                                          #
# miRNA associated with blood endophenotypes (excl. neutrophils)   #
#------------------------------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_6.png', width=750, height=750)
ggarrange(p_urea_covars, p_chloride_covars, p_crp_covars, 
          p_il6_covars, p_lymphocyte_covars, 
          p_nl_ratio_covars, p_ddimer_covars, 
          ncol=3, nrow = 3,
          labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
          label.x = 0, label.y = 1,
          font.label = list(size=18, color='black', face='bold'))
dev.off()


#--------------------------------------------#
# SUPPLEMENTARY TABLE 6: miRNAs associated   #
# with blood phenotypes                      #
#--------------------------------------------#
sig_table = rbind(sig_results_covars_df_urea, sig_results_covars_df_chloride, 
                  sig_results_covars_df_crp, sig_results_covars_df_il6,
                  sig_results_covars_df_neutrophil, sig_results_covars_df_lymphocyte,
                  sig_results_covars_df_nl_ratio,
                  sig_results_covars_df_ddimer)

sig_table$effect = round(as.numeric(sig_table$effect), 6)
sig_table = add_sig_stars(sig_table, 'p')

names(sig_table) = c('Dependent variable', 'miRNA', 'Covars', 'Effect', 
                     'P-value', 'N', 'R2 adjusted', 'Top SNP', 'Significance')

sig_table$miRNA = str_replace_all(sig_table$miRNA,
                                  '\\.', "-")

write.csv(sig_table, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_6A.csv')

sig_table_01 = sig_table %>% filter(`P-value` < 0.01)
n_unique_mirna = length(unique(sig_table_01$miRNA))

print(str_interp("There are ${n_unique_mirna} unique miRNAs associated with at least one blood pheno."))


#--------------------------------------------------------#
# SUPPLEMENTARY FIGURE 7:                                #
# Replication of miRNA associations for blood biomarkers #
#--------------------------------------------------------#
p_rep_urea = make_replication_plot(sig_results_df_urea,
                                   sig_results_covars_df_urea, 'Urea')
p_rep_cl = make_replication_plot(sig_results_df_chloride,
                                 sig_results_covars_df_chloride, 'Chloride')
p_rep_crp = make_replication_plot(sig_results_df_crp,
                                  sig_results_covars_df_crp, 'CRP')
p_rep_il6 = make_replication_plot(sig_results_df_il6,
                                  sig_results_covars_df_il6, 'IL6')
p_rep_neutro = make_replication_plot(sig_results_df_neutrophil,
                                     sig_results_covars_df_neutrophil, 'Neutrophil count')
p_rep_lympho = make_replication_plot(sig_results_df_lymphocyte,
                                     sig_results_covars_df_lymphocyte, 'Lymphocyte count')
p_rep_nl = make_replication_plot(sig_results_df_nl_ratio,
                                 sig_results_covars_df_nl_ratio, 'Neutrophil-to-lymphocyte ratio')
p_rep_ddimer = make_replication_plot(sig_results_df_ddimer,
                                     sig_results_covars_df_ddimer, 'D-dimers')


png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_7.png', width=900, height=900)
ggarrange(p_rep_urea, p_rep_cl,
          p_rep_crp, p_rep_il6, p_rep_neutro, p_rep_lympho, p_rep_nl,
          p_rep_ddimer,
          nrow=3, ncol=3, 
          labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'),
          label.x = 0, label.y = 1, 
          font.label = list(size=18, color='black', face='bold'))
dev.off()






