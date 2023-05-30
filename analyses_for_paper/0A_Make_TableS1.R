rm(list=ls())

library(dplyr)

#-----------#
# LOAD DATA # 
#-----------#
data = read.csv('~/Desktop/COVID_miRNA/Processed_data/data_merged_w_viral.csv')
variables = read.csv('~/Desktop/COVID_miRNA/Processed_data/Seha_data_variables.csv')

data_TableS1 = data %>%
  select(Barcode, ICU, Age_numeric, GENDER, onset_to_admit, 
         Region, BMI,
         FIRST_TEMPERATURE_ORAL, FIRST_OXYGEN_SATURATION, FIRST_RESPIRATORY_RATE,
         FIRST_SYSTOLIC, FIRST_DIASTOLIC, BLOOD_GROUP, 
         HYPERTENSION_YN, DIABETES_YN, ASTHMA_YN, ISCHEMIC_HEART_DISEASE_YN,
         CANCER_YN, ACUTE_KIDNEY_FAILURE_YN, CKD_YN, FUNGAL_INFECTION_YN, HIV_YN, 
         SEPSIS_YN, MYOCARDIAL_INFARCTION_YN, N_conditions, 
         Fever, Fatigue, Headache, Cough, SoreThroat, DifficultyBreathing, RunnyNose,
         MusclePain, Nausea, LossApetite, LossTaste, LossSmell, EyePain,
         FIRST_UREA_LVL, FIRST_CHLORIDE_LVL, FIRST_CO2, FIRST_POTASSIUM_LVL, FIRST_MAGNESIUM_LVL,
         FIRST_PHOSPHOROUS_LVL, FIRST_BILI_TOTAL, FIRST_C_REACTIVE_PROT, FIRST_TROPONIN_T,
         FIRST_VITAMIN_D_OH_LVL, FIRST_INTERLEUKIN6, FIRST_ABSOLUTE_NEUTROPHIL, FIRST_ABSOLUTE_LYMPHOCYTE,
         FIRST_NEUTRO_TO_LYMPHO_RATIO, FIRST_HBA1C, FIRST_PROCALCITONIN, FIRST_D_DIMER, FIRST_HEPB,
         log2_viral, miRNA)

write.csv(data_TableS1, '/Users/tamigjorgjieva/Desktop/covid_paper_outputs_oasis/Supplementary_Table_1.csv', 
          row.names=FALSE)
