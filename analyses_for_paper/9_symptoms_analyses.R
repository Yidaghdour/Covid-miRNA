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

# Functions for categorical variables
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

plot_volcano_logit = function(results_df, sig_h1, sig_h2, sig_l, dep_var, save, label=TRUE){
  
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
  
  
  
  if (save != 'no'){
    png(save, width=500, height=400)
    plot(p)
    dev.off()
  }
  
  return(p)
}


#------------#
# LOAD DATA  #
#------------#
setwd('~/Desktop/miRNA/mirna_processed/n96/')
data = read.csv('mirna_oasis_5_in_50_percent_standard_clinical.csv')

mirna_list = names(data)[3:634]

covars = c('Age_numeric', 'Gender_binary')

## Tables
table_outcomes = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4A.csv')
table_biomarkers = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_6A.csv')


#============================#
# 1. PHYSICAL MEASUREMENTS   #
#============================#

## TEMPERATURE -----------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_temp')

# Prepare dependent variable
data$FIRST_TEMPERATURE_ORAL[data$FIRST_TEMPERATURE_ORAL == 0] = NA
data$temp_scaled = scale(data$FIRST_TEMPERATURE_ORAL, center=TRUE, scale=TRUE)

# Compute OLS results 
results_covars_df = run_ols_covars('temp_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_covars_df_temp = return_and_save_table(results_covars_df, 0.05, 'table_temp_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_temp_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Body temperature', 'no', TRUE)


## OXYGEN SATURATION -----------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_oxy')

# Prepare dependent variable
data$FIRST_OXYGEN_SATURATION[data$FIRST_OXYGEN_SATURATION == 0] = NA
data$oxy_scaled = scale(data$FIRST_OXYGEN_SATURATION, center=TRUE, scale=TRUE)

# Compute OLS results 
results_covars_df = run_ols_covars('oxy_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_covars_df_oxy = return_and_save_table(results_covars_df, 0.05, 
                                                   'table_oxy_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_oxy_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Oxygen saturation', 'no', TRUE)


## RESPIRATORY RATE -----------------------------------------------------------------------------------
setwd('~/Desktop/mirna/gwas_oasis/mirna_resp')

# Prepare dependent variable
data$FIRST_RESPIRATORY_RATE[data$FIRST_RESPIRATORY_RATE == 0] = NA
data$resp_scaled = scale(data$FIRST_RESPIRATORY_RATE, center=TRUE, scale=TRUE)

# Compute OLS results 
results_covars_df = run_ols_covars('resp_scaled', mirna_list, TRUE, 3, covars)
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")

# Save results as a table 
sig_results_covars_df_resp = return_and_save_table(results_covars_df, 0.05, 'table_resp_scaled')

# Plot results (p<0.05 blue, p<0.01 red, top miRNA labeled)
p_resp_covars = plot_volcano(results_covars_df, 0.05, 0.01, 'Respiratory rate', 'no', TRUE)


#=========================#
# 2. NUMBER OF SYMPTOMS   #
#=========================#

## Count the number of self-reported symptoms at admission 
data_n_symptoms = data %>%
  select(Barcode, Fever:None)

data_n_symptoms$n_symptoms = rowSums(data_n_symptoms[,2:14])
mean = mean(na.omit(data_n_symptoms$n_symptoms))

# Median of number of symptoms = 4; Mean = 4.08

## Define highly symptomatic as having >3 self-reported symptoms
data_n_symptoms$highly_symptomatic = NA 

data_n_symptoms[is.na(data_n_symptoms$n_symptoms) == FALSE & 
                  data_n_symptoms$n_symptoms < mean, 'highly_symptomatic'] = 0
data_n_symptoms[is.na(data_n_symptoms$n_symptoms) == FALSE & 
                  data_n_symptoms$n_symptoms > mean, 'highly_symptomatic'] = 1

data = merge(data, data_n_symptoms %>% 
               select(Barcode, n_symptoms, highly_symptomatic), by='Barcode')

## Compute miRNAs associated with being highly symptomatic 

# Compute OLS logit results 
results_covars_df = run_ols_logit_covars('highly_symptomatic', mirna_list, covars)

# Save results as a table 
sig_results_covars_df_symptoms = return_and_save_table_logit(results_covars_df, 0.05, 'table_symptoms')

# Fix up table for plotting
results_covars_df$x = str_replace_all(results_covars_df$x, "\\.", "-")
results_covars_df$x = str_replace_all(results_covars_df$x, "hsa-", "")

# Plot results (p<0.05 highlighted, p<0.005 labeled)
p_symptoms_covar = plot_volcano_logit(results_covars_df, 0.05, 0.01, 0.005,
                                      'Highly symptomatic', 'no')


#-------------------------------------------#
# SUPPLEMENTARY TABLE 10A:                  #
# Associations between miRNA and symptoms   #
#-------------------------------------------#
sig_table = rbind(sig_results_covars_df_temp, sig_results_covars_df_resp,
                  sig_results_covars_df_oxy)

sig_table$effect = round(as.numeric(sig_table$effect), 6)
sig_table = add_sig_stars(sig_table, 'p')

names(sig_table) = c('Dependent variable', 'miRNA', 'Covars', 'Effect', 
                     'P-value', 'N', 'R2 adjusted', 'Top SNP', 'Significance')

sig_table$miRNA = str_replace_all(sig_table$miRNA, '\\.', "-")

write.csv(sig_table, "~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_10A.csv",
          row.names=FALSE)


#---------------------------------------------------------#
# SUPPLEMENTARY TABLE 10B:                                #
# Associations between miRNA and symptomatology category  #
#---------------------------------------------------------#
sig_results_covars_df_symptoms = add_sig_stars(sig_results_covars_df_symptoms, 'p')
names(sig_results_covars_df_symptoms) = c('Dependent variable', 'miRNA', 'Covars', 'Effect', 
                                          'P-value', 'N', 'AIC', 'Significance')

sig_results_covars_df_symptoms$miRNA = str_replace_all(sig_results_covars_df_symptoms$miRNA,
                                                       '\\.', "-")

write.csv(sig_results_covars_df_symptoms, 
          "~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_10B.csv",
          row.names=FALSE)


#-----------------------------------------------------#
# SUPPLEMENTARY FIGURE 11:                            #
# miRNAs associated with severity of symptomatology   #
#-----------------------------------------------------#
p_blank = ggplot() + theme_void()

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_11.png', width=750, height=550)
ggarrange(p_temp_covars, p_resp_covars, p_oxy_covars,
          p_blank, p_symptoms_covar, p_blank,
          labels = c('A', 'B', 'C', '', 'D', ''),
          ncol=3, nrow = 2,
          label.x = 0, label.y = 1,
          font.label = list(size=18, color='black', face='bold'))
dev.off()


