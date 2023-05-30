rm(list=ls())

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(MetBrewer)

#-----------#
# FUNCTIONS # 
#-----------#

make_freq_table = function(data, col_name, order=NA){
  
  freq_df = data.frame(table(data[,col_name]))

  # Order by provided list or otherwise, by frequency
  if (sum(is.na(order)) != 0){
    freq_df = freq_df[order(-freq_df$Freq),]
  } else {
    freq_df = freq_df[order(match(freq_df$Var1, order)), ]
  }
  
  freq_df$Percentage = paste(
    round(freq_df$Freq / nrow(data)*100, 1), "%", sep='')
  
  return(freq_df)
  
}

make_freq_based_on_multiple_cols = function(colnames_list, order, yes, min_freq=5, group){
  
  # Make empty table
  data_summary = data.frame(matrix(nrow=0, ncol=5))
  names(data_summary) = c('Characteristic', 'Freq (full)', 'Percent (full)',
                          'Freq (miRNA sample)', 'Percent (miRNA sample)')
  
  for (var in colnames_list){
    
    # Make frequency tables for full dataset (N=259) and miRNA data (N=96)
    freq_table_full = make_freq_table(data, var, order)
    freq_table_mirna = make_freq_table(data_mirna, var, order)
    
    # Merge the two tables 
    freq_tables = merge(freq_table_full, freq_table_mirna,
                        by='Var1', all=TRUE, sort=FALSE)
    names(freq_tables) = names(data_summary)
    
    # Keep only 'Yes' and rename 'Characteristic' to var
    if (nrow(freq_tables) == 2){
      freq_tables = freq_tables[freq_tables$Characteristic == yes,]
      
      if (group == 'Conditions'){
        freq_tables$Characteristic = paste(substr(var, 1, 1), 
                                           tolower(substr(var, 2, nchar(var))), sep='')
      } else {
        freq_tables$Characteristic = var
      }
        
      # Add to master table if frequency > 5
      if (freq_tables$`Freq (full)` > min_freq){
        data_summary = rbind(data_summary, freq_tables)
      }
    }
    
  }
  
  data_summary = data_summary[order(-data_summary$`Freq (full)`),]
  
  return(data_summary)
}

replace_values = function(df, col_name, replacements_df){
  
  df[,col_name] = as.character(df[,col_name])
  
  # Replace values from 1st column of replacements_df with values in 2nd column
  for (i in 1:nrow(replacements_df)){
    
    to_replace = replacements_df[i,1]
    replace_with = replacements_df[i,2]
    
    df[,col_name][df[,col_name] == to_replace] = replace_with
  }
  
  return(df)
  
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



#-----------#
# LOAD DATA # 
#-----------#
data = read.csv('~/Desktop/COVID_miRNA/Processed_data/data_merged_w_viral.csv')
variables = read.csv('~/Desktop/COVID_miRNA/Processed_data/Seha_data_variables.csv')

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

# ICU
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


#===============================================================#
# TABLE 1:                                                      #
# Demographics (gender, age, region) for full dataset (n=259)   #
# and the miRNA dataset (n=96)                                  #
#===============================================================#

data_summary = data.frame(matrix(nrow=0, ncol=5))
names(data_summary) = c('Characteristic', 'Freq (full)', 'Percent (full)',
                        'Freq (miRNA sample)', 'Percent (miRNA sample)')

vars = c('GENDER', 'Age_bins', 'Region')
orders = list(gender_factors, age_factors, regions_factors)

for (i in 1:length(vars)){
  
  # Select variable and order of variables in table 
  var = vars[i]
  order = orders[[i]]

  # Make frequency tables for full dataset (N=259) and miRNA data (N=96)
  freq_table_full = make_freq_table(data, var, order)
  freq_table_mirna = make_freq_table(data_mirna, var, order)
  
  # Merge the two tables 
  freq_tables = merge(freq_table_full, freq_table_mirna,
                      by='Var1', all=TRUE, sort=FALSE)
  names(freq_tables) = names(data_summary)
  
  # Add to master table 
  data_summary = rbind(data_summary, freq_tables)
}

data_summary = replace_values(data_summary, 'Characteristic', regions_df)
write.csv(data_summary, '~/Desktop/covid_paper_outputs_oasis/Table_1.csv')


#=====================================================#
# SUPPLEMENTARY TABLE 1: Pre-existing conditions      #
# for full dataset (n=259) and miRNA dataset (n=96)   #
#=====================================================#

# Identify comorbidities variables 
conditions_vars = variables[variables$Category == 'Medical conditions', 
                            'Variable']
conditions_vars = conditions_vars[!(endsWith(conditions_vars, '_YN'))]
conditions_vars = conditions_vars[conditions_vars != "N_conditions"]

# Make frequency table
data_summary = make_freq_based_on_multiple_cols(conditions_vars, icu_factors, 'Yes', 5,
                                                'Conditions')

write.csv(data_summary, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_1.csv')


#=====================================================#
# SUPPLEMENTARY TABLE 2: Symptoms at admission for    #
# full dataset (n=259) and miRNA dataset (n=96)       #
#=====================================================#

# Identify symptom variables
symptoms_vars = names(data %>%
  dplyr::select(Fever:EyePain))

# Make frequency table
data_summary = make_freq_based_on_multiple_cols(symptoms_vars, NA, 1, 5, 'Not')

# Replace symptoms with better text 
data_summary = replace_values(data_summary, 'Characteristic', symptoms_df)

write.csv(data_summary, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_2.csv')
