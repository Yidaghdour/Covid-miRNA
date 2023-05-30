rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(ggrepel)


#-----------#
# FUNCTIONS # 
#-----------#
make_var_effect_and_sig = function(tables){
  tables$Effect_and_sig = 0
  
  tables[tables$Effect < 0 & tables$Significance == '****', 'Effect_and_sig'] = 1
  tables[tables$Effect < 0 & tables$Significance == '***', 'Effect_and_sig'] = 2
  tables[tables$Effect < 0 & tables$Significance == '**', 'Effect_and_sig'] = 3
  tables[tables$Effect < 0 & tables$Significance == '*', 'Effect_and_sig'] = 4
  
  tables[tables$Effect > 0 & tables$Significance == '*', 'Effect_and_sig'] = 5
  tables[tables$Effect > 0 & tables$Significance == '**', 'Effect_and_sig'] = 6
  tables[tables$Effect > 0 & tables$Significance == '***', 'Effect_and_sig'] = 7
  tables[tables$Effect > 0 & tables$Significance == '****', 'Effect_and_sig'] = 8
  
  return(tables)
}

factorize_effect_and_sig = function(tables, ann_table){
  
  # Factorize the effect + significance
  tables$Effect_and_sig_factor = tables$Effect_and_sig
  tables$Effect_and_sig_factor = factor(tables$Effect_and_sig)
  
  levels_num = levels(tables$Effect_and_sig_factor)
  levels_txt = c()
  
  for (level_num in levels_num){
    level_txt = ann_table[which(ann_table$Numeric == level_num), 'Text']
    levels_txt = c(levels_txt, level_txt)
  }
  
  levels(tables$Effect_and_sig_factor) = levels_txt
  
  return(tables)
}


#--------#
# LISTS  # 
#--------#
var_names = c("out_ICU", "urea_scaled", "chloride_scaled", 
              "crp_scaled", "il6_scaled", "neutrophil_scaled", "lymphocyte_scaled",
              'neutro_lympho_ratio_scaled',
              "ddimer_scaled")

var_names_txt = c('ICU admission', 'Urea', 'Chloride', 
                  'CRP', 'IL6', 'Neutrophil count', 'Lymphocyte count', 'Neutro/lympho ratio',
                  'D-dimers')

var_df = data.frame(cbind(var_names, var_names_txt))
names(var_df) = c('Variable names', 'Variable text')

ann_table = data.frame(cbind(c(1:8),
                             c('Negative effect, p < 0.0001', 'Negative effect, p < 0.001',
                               'Negative effect, p < 0.01', 'Negative effect, p < 0.05', 
                               'Positive effect, p < 0.05', 'Positive effect, p < 0.01',
                               'Positive effect, p < 0.001', 'Positive effect, p < 0.0001')))
names(ann_table) = c('Numeric', 'Text')


#-----------#
# LOAD DATA # 
#-----------#
table_outcomes = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4A.csv')
table_biomarkers = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_6.csv')


#---------------------------#
# DATA PREP & VISUALIZATION # 
#---------------------------#

# Prepare tables for merging 
table_outcomes = table_outcomes %>% 
  dplyr::select(Dependent.variable, miRNA, Effect, P.value, Significance)

table_biomarkers = table_biomarkers %>% 
  dplyr::select(Dependent.variable, miRNA, Effect, P.value, Significance)

# Merge all results in one table
tables = rbind(table_outcomes, table_biomarkers)

# Make a variable for effect AND significance
tables = make_var_effect_and_sig(tables)

# Factorize the effect and significance column
tables = factorize_effect_and_sig(tables, ann_table)

# Rename dependent variables
tables$Dependent.variable = factor(tables$Dependent.variable,
                                   levels = var_names)
levels(tables$Dependent.variable) = var_names_txt

tables = tables %>%
  filter(P.value < 0.01)

#------------------------------------------------------------------#
# SUPPLEMENTARY FIGURE 7:                                          #
# miRNAs associated with ICU are associated with blood phenotypes  #
#------------------------------------------------------------------#

# Prepare panel 
icu_mirna = table_outcomes[table_outcomes$Dependent.variable == 'out_ICU' &
                             table_outcomes$P.value < 0.01, 'miRNA']

tables_shared_mirna = tables %>%
  filter(miRNA %in% icu_mirna)

tmp_df = tables_shared_mirna %>%
  filter(Dependent.variable == 'ICU admission')

mirna_ordered = tmp_df[order(tmp_df$Effect_and_sig), 'miRNA']
tables_shared_mirna$miRNA = factor(tables_shared_mirna$miRNA,
                                   levels = mirna_ordered)

tmp_df = tables_shared_mirna %>% 
  filter(Dependent.variable != 'ICU admission')

n_mirnas_icu_and_blood_pheno = length(unique(tmp_df$miRNA))
n_mirnas_icu = length(icu_mirna)
perc = round((n_mirnas_icu_and_blood_pheno / n_mirnas_icu)*100, 1)

print(str_interp("${n_mirnas_icu_and_blood_pheno} out of ${n_mirnas_icu} (${perc}%) miRNAs associated with ICU are also associated with blood phenotypes."))

# Prepare colors
factors_i = sort(unique(tables$Effect_and_sig))
colors = rev(brewer.pal(n = 8, name = "RdBu"))[factors_i]

# Plot
p = ggplot(tables_shared_mirna, aes(x=factor(Dependent.variable), y=(factor(miRNA)), 
                                      fill=factor(Effect_and_sig_factor))) +
  geom_tile(color='black', size=0.3, width=1) +
  theme_classic() +
  theme(axis.text.x = element_text(size=13, angle=60, vjust=1, hjust=1),
        axis.text.y = element_text(size=13),
        legend.title = element_blank(),
        plot.title = element_text(size=15, face='bold', hjust=0.5),
        legend.text = element_text(size=12),
        plot.margin = margin(0.5, 0.5, 0, 0, 'cm'),
        legend.position = 'bottom') +
  ylab('') + xlab('') +
  scale_fill_manual(values = colors) +
  ggtitle("miRNAs associated with ICU admission") +
  guides(fill=guide_legend(nrow=4, byrow=TRUE))

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_7.png', width=500, height=700)
plot(p)
dev.off()


#--------------------------------------------------#
# SUPPLEMENTARY FIGURE 11:                         #
# Shared miRNA associations across endophenotypes  #
#--------------------------------------------------#

# Identify miRNAs associated with multiple phenotypes 
freq_df = data.frame(table(tables$miRNA))

freq_df = freq_df[freq_df$Freq > 3, ]
shared_mirna_ordered = freq_df[order(freq_df$Freq), 'Var1']

# Subset from tables and visualize
shared_mirna_df = tables %>%
  filter(miRNA %in% shared_mirna_ordered)

shared_mirna_df$miRNA = factor(shared_mirna_df$miRNA, 
                              levels = shared_mirna_ordered)

factors_i = sort(unique(shared_mirna_df$Effect_and_sig))
colors = rev(brewer.pal(n = 8, name = "RdBu"))[factors_i]

# Plot
p = ggplot(shared_mirna_df, aes(x=factor(Dependent.variable), y=factor(miRNA), 
                           fill=factor(Effect_and_sig_factor))) +
  geom_tile(color='black', size=0.3, width=1) +
  theme_classic() +
  geom_text(aes(label = Significance), size=5.5) +
  theme(axis.text.x = element_text(size=14, angle=60, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=12)) +
  ylab('') + xlab('') +
  guides(fill=guide_legend(nrow=4, byrow=TRUE)) +
  scale_fill_manual(values = colors)

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_11.png', width=450, height=600)
plot(p)
dev.off()


#----------------------------------------------#
# SUPPLEMENTARY FIGURE 12:                     #
# Number of shared miRNAs across all factors   #
#----------------------------------------------#

factors = unique(as.character(tables$Dependent.variable))

summary_df = data.frame(matrix(nrow=0, ncol=4))
names(summary_df) = c('Factor1', 'Factor2', 'miRNAs in common', 'N_mirna')

for (i in 1:(length(factors)-1)){

  factor1 = factors[i]
  other_factors = factors[(i+1):length(factors)]
  
  for (factor2 in other_factors){

    mirna_factor1 = tables[tables$Dependent.variable == factor1, 'miRNA']
    mirna_factor2 = tables[tables$Dependent.variable == factor2, 'miRNA']
    
    mirna_in_common = intersect(mirna_factor1, mirna_factor2)
    mirna_in_common_txt = paste(mirna_in_common, collapse=', ')
    mirna_in_common_n = length(mirna_in_common)
    
    summary_df[nrow(summary_df)+1, ] = c(factor1, factor2, mirna_in_common_txt, mirna_in_common_n)
  }
}

summary_df[summary_df == 0] = NA

summary_df$Factor1 = factor(summary_df$Factor1, levels = (factors))
summary_df$Factor2 = factor(summary_df$Factor2, levels = rev(factors))
summary_df$N_mirna = as.numeric(summary_df$N_mirna)

p = ggplot(summary_df, aes(x=Factor1, y=Factor2, fill=N_mirna)) +
  geom_tile() +
  geom_text(aes(label=N_mirna), size=5) +
  theme_classic() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1),
        legend.position = 'none') +
  xlab('') + ylab('') + 
  scale_fill_distiller(palette='Reds', na.value = 'white', trans = "reverse")

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_12.png', width=430, height=430)
plot(p)
dev.off()


#---------------------------------------#
# SUPPLEMENTARY FIGURE 13:              #
# Pairs of factors with shared miRNAs   #
#---------------------------------------#
all_factors_wide = tables %>% 
  select(miRNA, Dependent.variable, Effect) %>%
  spread(Dependent.variable, Effect)

summary_df_subset = summary_df %>%
  filter(N_mirna >= 5) %>%
  filter(Factor1 != 'ICU admission') %>%
  filter(Factor2 != 'ICU admission') %>%
  filter(Factor1 != 'Oxygen therapy') %>%
  filter(Factor2 != 'Oxygen therapy')

count=1

for (i in 1:nrow(summary_df_subset)){

  f1_txt = as.character(summary_df_subset[i,'Factor1'])
  f2_txt = as.character(summary_df_subset[i,'Factor2'])

  tmp = all_factors_wide[,c(f1_txt, f2_txt)]
  names(tmp) = c('f1', 'f2')
  
  min = min(min(na.omit(tmp[,1])), min(na.omit(tmp[,2])))
  max = max(max(na.omit(tmp[,1])), max(na.omit(tmp[,2])))
  
  total_max = max(abs(min), abs(max))
  total_tick = floor(total_max)
  
  min = min(min(na.omit(tmp[,1])), min(na.omit(tmp[,2])))
  max = max(max(na.omit(tmp[,1])), max(na.omit(tmp[,2])))
  
  total_max = max(abs(min), abs(max))
  
  assign(str_interp("p${count}"),
         ggplot(tmp, aes(x=f1, y=f2)) +
           geom_point(shape=21, color='black', fill='#3FA796', size=4) +
           theme_classic() +
           theme(axis.text = element_text(size=14),
                 axis.title = element_text(size=14),
                 plot.title = element_text(size=14, face='bold', hjust=0.5)) +
           xlab('') +
           ylab('') +
           geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
           geom_hline(yintercept=0, linetype='dashed', color='grey') +
           geom_vline(xintercept=0, linetype='dashed', color='grey') +
           ggtitle(str_interp("${f1_txt} vs\n${f2_txt}")) +
           scale_x_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick)) +
           scale_y_continuous(limits = c(-total_max, total_max), breaks=c(-total_tick, 0, total_tick)))
  
  count = count + 1
}

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_13.png', width=600, height=800)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7,
                        p8, p9, p10, nrow=4)
dev.off()

