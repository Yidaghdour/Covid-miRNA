rm(list=ls())

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(MetBrewer)
library(corrplot)


#-----------#
# FUNCTIONS # 
#-----------#
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

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

add_mt_correction_pearson_p = function(df, correction_method, adj_p_colname){
  
  df$cor_pearson_p = as.numeric(df$cor_pearson_p)
  
  df$Adjusted_P = p.adjust(df[,'cor_pearson_p'], 
                      method=correction_method, n=nrow(df))
  
  n_sig_nominal = nrow(df %>% filter(cor_pearson_p < 0.05))
  n_sig_adjusted = nrow(df %>% filter(Adjusted_P < 0.05))
  
  names(df)[ncol(df)] = adj_p_colname
  
  print(str_interp("There are ${n_sig_nominal} associations with nominal P < 0.05."))
  print(str_interp("There are ${n_sig_adjusted} associations with ${correction_method} P < 0.05."))
  
  return(df)
  
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


#-----------#
# LOAD DATA # 
#-----------#
data = read.csv('~/Desktop/COVID_miRNA/Processed_data/data_merged_w_viral.csv')
variables = read.csv('~/Desktop/COVID_miRNA/Processed_data/Seha_data_variables.csv')

# Make neutrophil-to-lymphocyte ratio variable
data$FIRST_NEUTRO_TO_LYMPHO_RATIO = data$FIRST_ABSOLUTE_NEUTROPHIL/data$FIRST_ABSOLUTE_LYMPHOCYTE
data$FIRST_NEUTRO_TO_LYMPHO_RATIO[data$FIRST_NEUTRO_TO_LYMPHO_RATIO %in% c('NaN', 'Inf')] = NA

# Add this to variables
variables = rbind(variables[1:336, ], 
                  c('FIRST_NEUTRO_TO_LYMPHO_RATIO', 'Lab tests', 'numeric', 'Continuous', 'Not'),
                  variables[337:nrow(variables), ])

# Subset miRNA dataset
data_mirna = data %>% 
  filter(miRNA == 1)


#-----------#
# LISTS     # 
#-----------#

# Age
age_factors = c('18-25 yrs', '25-34 yrs', '35-44 yrs', '45-54 yrs',
                '55-64 yrs', '65-74 yrs', '75-84 yrs', '85+ yrs')

# Gender
gender_factors = c('Male', 'Female')

# Regions
regions_factors = c('MENA', 'Africa_excl_MENA', 'Southeast_Asia', 
                    'NorthAmerica_Australia')
regions_txt = c('MENA (Middle East and North Africa',
                'Africa (excluding MENA)',
                'Southeast Asia',
                'North America and Australia')
regions_df = data.frame(cbind(regions_factors, regions_txt))
names(regions_df) = c('Factor', 'Text')

# BMI categories
bmi_factors = c('Underweight', 'Normal', 'Obese', 'Overweight')
bmi_txt = c('Underweight (BMI < 18.5)', 
            'Normal or healthy weight (18.5 - 24.9)',
            'Overweight (25 - 29.9)',
            'Obese (30 and above)')

bmi_df = data.frame(cbind(bmi_factors, bmi_txt))
names(bmi_df) = c('Factor', 'Text')

# Blood groups categories
blood_factors = c('O POS', 'O NEG', 'A POS', 'A NEG',
                  'B POS', 'B NEG', 'AB POS', 'AB NEG')
blood_txt = c('O positive', 'O negative', 'A positive', 'A negative',
              'B positive', 'B negative', 'AB positive', 'AB negative')
blood_df = data.frame(cbind(blood_factors, blood_txt))
names(blood_df) = c('Factor', 'Text')

# ICU factors
icu_factors = c('Yes', 'No')

# Patient days
patient_days_factors = c('0-5 days', '5-10 days', '10-20 days', '20-30 days',
                         '30-40 days', '40-50 days', '50+ days')

# Symptoms
symptoms_factors = names(data %>%
                           dplyr::select(Fever:EyePain))
symptoms_txt = c('Fever', 'Fatigue', 'Headache', 'Cough', 'Sore throat',
                 'Difficulty breathing', 'Runny nose', 'Muscle pain', 'Nausea',
                 'Loss of apetite', 'Loss of taste', 'Loss of smell', 'Eye pain')
symptoms_df = data.frame(cbind(symptoms_factors, symptoms_txt))
names(symptoms_df) = c('Factor', 'Text')

# Predictors for correlation and OLS analyses
demographic_vars = c('Age_numeric', 'Gender_binary', 'MENA_YN', 'Southeast_Asia_YN',
                     'Africa_excl_MENA_YN', 'NorthAmerica_Australia_YN')

factors_vars = c('BMI', 'blood_A', 'blood_B', 'blood_AB', 'blood_O', 'blood_Pos', 'blood_Neg')

conditions_vars = c('HYPERTENSION_YN', 'DIABETES_YN', 'ASTHMA_YN', 'ISCHEMIC_HEART_DISEASE_YN', 
                    'CANCER_YN', 'ACUTE_KIDNEY_FAILURE_YN', 'CKD_YN', 'FUNGAL_INFECTION_YN', 
                    'HIV_YN', 'SEPSIS_YN', 'MYOCARDIAL_INFARCTION_YN', 'N_conditions')

symptoms_vars = symptoms_factors

physical_measurements = c('FIRST_TEMPERATURE_ORAL', 'FIRST_OXYGEN_SATURATION', 
                          'FIRST_RESPIRATORY_RATE', 'FIRST_SYSTOLIC', 'FIRST_DIASTOLIC')

first_lab_tests = variables[variables$Category == 'Lab tests' & variables$Class == 'numeric' & 
                              startsWith(variables$Variable, 'FIRST_'), 'Variable']

all_predictors = c(demographic_vars, factors_vars,
                   conditions_vars, symptoms_vars,
                   physical_measurements, first_lab_tests, 'log2_viral')

all_predictors_txt = c('Age', 'Gender', 
                       'Middle East or North Africa', 
                       'Southeast Asia',
                       'Africa (excluding MENA)',
                       'North America and Australia',
                       'BMI',
                       'Blood group A', 'Blood group B',
                       'Blood group AB', 'Blood group O',
                       'Blood group - positive', 
                       'Blood group - negative',
                       'Hypertension', 'Diabetes', 'Asthma',
                       'Ischemic heart disease', 'Cancer',
                       'Acute kidney failure', 'Chronic kidney disease',
                       'Fungal infection', 'HIV', 'Sepsis',
                       'Myocardial infarction', 'Number of existing conditions',
                       'Fever', 'Fatigue', 'Headache', 'Cough', 
                       'Sore throat', 'Difficulty breathing', 'Runny nose', 
                       'Muscle pain', 'Nausea', 'Loss of apetite', 
                       'Loss of sense of taste', 'Loss of sense of smell', 'Eye pain',
                       'Temperature', 'Oxygen saturation', 
                       'Respiratory rate', 'Systolic pressure',
                       'Diastolic pressure', 'Urea', 'Chloride', 'CO2',
                       'Potassium', 'Magnesium', 'Phosphorous', 'Bilirubin',
                       'C reactive protein (CRP)',
                       'Troponin T', 'Vitamin D OH', 'Inteleukin 6 (IL6)',
                       'Absolute neutrophil number',
                       'Absolute lymphocyte number', 'Neutrophil-to-lymphocyte ratio',
                       'HbA1c', 'Procalcitonin', 'D-dimers', 'HEPB', 
                       'Viral load')

categories = c(rep('Demographics', length(demographic_vars)),
               rep('Factors', length(factors_vars)),
               rep('Conditions', length(conditions_vars)),
               rep('Symptoms', length(symptoms_vars)),
               rep('Physical measurements', length(physical_measurements)),
               rep('First lab tests', length(first_lab_tests)),
               'Viral load')


#--------------------------#
# ASSOCIATION ANALYSES     # 
#--------------------------#
predictors_df = data.frame(cbind(c(all_predictors),
                                 c(all_predictors_txt),
                                 c(categories),
                                 c(rep('out_ICU', length(all_predictors))),
                                 c(rep('ICU admission', length(all_predictors)))))

names(predictors_df) = c('Variable', 'Variable_txt', 'Category', 'Outcome', 'Outcome_txt')

predictors_df = predictors_df %>% add_column(
  cor_pearson = NA, cor_pearson_p = NA,
  cor_spearman = NA, cor_spearman_p = NA,
  cor_kendal = NA, cor_kendal_p = NA,
  ols_beta = NA, ols_se = NA, ols_r2 = NA, ols_p = NA, 
  ols_age_sex_beta = NA, ols_age_sex_se = NA, ols_age_sex_r2 = NA, ols_age_sex_p = NA,
  ols_age_sex_time_beta = NA, ols_age_sex_time_se = NA, 
  ols_age_sex_time_r2 = NA, ols_age_sex_time_p = NA)


for (i in 1:nrow(predictors_df)){
  
  df_list = c()
  
  # Calculate correlations ==========================================================
  predictor = predictors_df[i, 'Variable']
  outcome = predictors_df[i, 'Outcome']
  
  # Replace 0's with NA =============================================================
  if (startsWith(predictor, 'FIRST_')){
    data[,predictor][data[,predictor] == 0] = NA
  }
  
  for (method in c('pearson', 'spearman', 'kendal')){
    cor = cor.test(data[,predictor], data[,outcome], method=method, use='pairwise.complete.obs')
    df_list = c(df_list, round(as.numeric(cor$estimate[1]), 3), cor$p.value)
  }
  
  # Calculate OLS effects  =========================================================
  formula = str_interp("${outcome} ~ ${predictor}")
  formula_age_sex = str_interp("${outcome} ~ ${predictor} + Age_numeric + Gender_binary")
  formula_age_sex_time = str_interp("${outcome} ~ ${predictor} + Age_numeric + Gender_binary +
                                      onset_to_admit")
  
  results = summary(lm(formula, data=data))
  df_list = c(df_list, round(c(results$coefficients[2,1],
                               results$coefficients[2,2],
                               results$adj.r.squared,
                               results$coefficients[2,4]), 6))
  
  if (!(predictor %in% c('Age_numeric', 'Gender_binary'))){
    results = summary(lm(formula_age_sex, data=data))
    df_list = c(df_list, round(c(results$coefficients[2,1],
                                 results$coefficients[2,2],
                                 results$adj.r.squared,
                                 results$coefficients[2,4]), 6))
    
    if (startsWith(predictor, 'FIRST_') | predictor %in% symptoms_factors){
      results = summary(lm(formula_age_sex_time, data=data))
      df_list = c(df_list, round(c(results$coefficients[2,1],
                                   results$coefficients[2,2],
                                   results$adj.r.squared,
                                   results$coefficients[2,4]), 6))
    } else {
      df_list = c(df_list, rep(NA, 4))
    }
    
  } else {
    df_list = c(df_list, rep(NA, 8))
  }
  
  i_col = which(names(predictors_df) == 'cor_pearson')
  predictors_df[i,i_col:ncol(predictors_df)] = df_list
  
}

predictors_df$Variable_txt = factor(predictors_df$Variable_txt,
                                    levels=all_predictors_txt)

# FDR adjustment 
predictors_df = add_mt_correction_pearson_p(predictors_df, 'fdr', 'pearson_p_FDR')
predictors_df = add_sig_stars(predictors_df, 'pearson_p_FDR', 'Significance_FDR')

#---------------------------------------------------------------#
# FIGURE 1B:                                                    #
# Factors associated with COVID outcomes (ICU, oxygen therapy)  #
#---------------------------------------------------------------#
tmp_sig = predictors_df %>%
  filter(pearson_p_FDR < 0.05)

tmp_sig$Outcome_txt = "ICU\nadmission"

p = ggplot(tmp_sig, aes(Outcome_txt, Variable_txt, fill=cor_pearson)) + 
  geom_tile(color='white', size=1.5) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic() +
  geom_text(aes(label = round(cor_pearson, 2)), size=6) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=16),
        legend.position = 'none') +
  ylab('') + xlab('') +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(position = "top")

png('~/Desktop/covid_paper_outputs_oasis/Figure_1B.png', width=350, height=700)
plot(p)
dev.off()

p = ggplot(tmp_sig, aes(Outcome_txt, Variable_txt, fill=cor_pearson)) + 
  geom_tile(color='white', size=1.5) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic() +
  geom_text(aes(label = round(cor_pearson, 2)), size=10) +
  theme(axis.text.x = element_text(size=26),
        axis.text.y = element_text(size=28),
        legend.position = 'none') +
  ylab('') + xlab('') +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(position = "top")

png('~/Desktop/covid_paper_outputs_oasis/Figure_1B_hq.png', width=550, height=980)
plot(p)
dev.off()


#------------------------------------------------------------------#
# SUPPLEMENTARY FIGURE 1:                                          #
# Correlation structure of factors associated with COVID outcomes  #
#------------------------------------------------------------------#
vars = c(unique(tmp_sig$Variable), 'out_ICU')

data_subset = data[,vars]
names(data_subset) = c('South Asia', 'Acute kidney failure',
                       'Sepsis', 'Myocardial infarction',
                       'Fever', 'Cough', 
                       'Temperature', 'Oxygen saturation', 'Respiratory rate',
                       'Urea', 'Chloride', 'CRP', 'IL-6', 'Neutrophil count',
                       'Lymphocyte count', 'Neutrophil-to-lymphocyte ratio', 
                       'D-dimer', 'ICU admission')

cor_m = cor(data_subset, use='pairwise.complete.obs')
cor_p = cor.mtest(data_subset)

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_1.png', width=750, height=750)
corrplot(cor_m, type="upper", order="hclust", method='color',
         tl.col="black", tl.srt=90, 
         diag=FALSE, addCoef.col = 'black',
         p.mat = cor_p, sig.level = 0.05, insig = "blank")
dev.off()


#--------------------------------------------#
# SUPPLEMENTARY TABLE 3A: Correlations and   #
# OLS results for ICU-associated factors     #
#--------------------------------------------#

icu_sig = predictors_df %>%
  filter(pearson_p_FDR < 0.05) %>% 
  filter(outcome == 'out_ICU') %>%
  dplyr::select(Variable_txt, Outcome_txt,
                cor_pearson, cor_pearson_p, pearson_p_FDR,
                ols_beta, ols_se, ols_p, 
                ols_r2, ols_age_sex_r2, ols_age_sex_time_r2)

icu_sig = add_sig_stars(icu_sig, 'pearson_p_FDR', 'Significance')
icu_sig = icu_sig[,c(1:5, 12, 9:11)]

names(icu_sig) = c('Independent variable', 'Dependent variable', 
                   'Pearson correlation', 'P-value',
                   'FDR-adjusted P-value', 'Significance', 
                   'R2 (no covariates)', 'R2 (covars: age and sex)', 
                   'R2 (covars: age, sex, and time from symptom onset to admission)')

icu_sig = icu_sig[,1:6]

write.csv(icu_sig, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_3A.csv')

#------------------------------------------------#
# SUPPLEMENTARY TABLE 3B:                        #
# Frequencies and mean values in ICU vs non-ICU  #
#------------------------------------------------#
vars = (predictors_df %>%
          filter(pearson_p_FDR < 0.05))$Variable

data_icu = data %>% filter(out_ICU == 1)
data_not_icu = data %>% filter(out_ICU != 1)

summary_df = data.frame(matrix(nrow=0, ncol=5))
names(summary_df) = c('Variable', 'ICU_mean', 'ICU_sd', 'nonICU_mean', 'nonICU_sd')

for (var in vars){
  n_unique = length(na.omit(unique(data[,var])))
  
  n_icu = length(na.omit(data_icu[,var]))
  n_nonicu = length(na.omit(data_not_icu[,var]))
  
  if (n_unique == 2){
    
    freq_icu = round((sum(na.omit(data_icu[,var]) == 1)/n_icu)*100, 1)
    freq_nonicu = round((sum(na.omit(data_not_icu[,var]) == 1)/n_nonicu)*100, 1)
    
    summary_df[nrow(summary_df)+1, ] = c(var, freq_icu, NA, freq_nonicu, NA)
    
  } else {
    mean_icu = round(mean(na.omit(data_icu[,var])), 1)
    mean_nonicu = round(mean(na.omit(data_not_icu[,var])), 1)
    
    sd_icu = round(sd(na.omit(data_icu[,var])), 1)
    sd_nonicu = round(sd(na.omit(data_not_icu[,var])), 1)

    summary_df[nrow(summary_df)+1, ] = c(var, mean_icu, sd_icu, mean_nonicu, sd_nonicu)
  }
  
}

summary_df = replace_names(summary_df, 'Variable', all_predictors, all_predictors_txt)
summary_df[is.na(summary_df)] = ''

write.csv(summary_df, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_3B.csv')



