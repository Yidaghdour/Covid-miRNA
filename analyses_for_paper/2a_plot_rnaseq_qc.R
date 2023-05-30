rm(list=ls())

library(tidyr)
library(stringr)
library(ggplot2)
library(tidyselect)
library(dplyr)
library(ggpubr)


#------------#
# FUNCTIONS  #
#------------#
return_and_save_density_plot = function(data, before_after, palette, norm_type=NULL, extra=NULL){
  
  indiv_ids = names(data)
  n_indiv = length(indiv_ids)
  
  data_long = data %>% 
    gather(indiv_ID, mirna_Expression, all_of(indiv_ids))
  
  n = length(indiv_ids)
  
  if (before_after == 'before'){
    plot_title = paste(paste("Before normalization"))
    file_name = paste("density_before_normalization", "_n_", n, ".png", sep='')
  } else {
    plot_title = paste(paste(norm_type, "normalization"))
    file_name = paste("density_", norm_type, "_normalization_", before_after, "_n_", n, 
                      ".png", sep='')
  }
  
  p = ggplot(data_long, aes(x=mirna_Expression, color=indiv_ID)) +
    geom_density() + theme(legend.position='none') + theme_light() +
    ggtitle(plot_title) + 
    theme(axis.title = element_text(size=13),
          axis.text = element_text(size=12),
          plot.title = element_text(size=14, hjust=0.5, face="bold"),
          legend.position='none') +
    xlab("\nTranscript expression") + ylab("Density\n") + 
    scale_color_manual(values=get_palette(palette = palette, n_indiv))
  
  png(file_name, width=450, height=320)
  plot(p)
  dev.off()
  
  return(p)
}


#------------#
# LOAD DATA  #
#------------#
data_transformed = read.csv('~/Desktop/miRNA/rnaseq/rnaseq_data_transformed_smallest.csv')
data_normalized = read.csv('~/Desktop/miRNA/rnaseq/rnaseq_data_normalized_smallest.csv')

rownames(data_transformed) = data_transformed$X
rownames(data_normalized) = data_normalized$X

data_transformed = data_transformed[,-1]
data_normalized = data_normalized[,-1]


#----------------#
# MAKE QC PLOTS  #
#----------------#
p_before = return_and_save_density_plot(data_transformed, 
                                        'before', 'Oranges')
p_after = return_and_save_density_plot(data_normalized, 
                                       'after', 'Purples', 'Standard')


#---------------------------------------------------#
# SUPPLEMENTARY FIGURE 4:                           #
# Before and after normalization for RNA-seq data   #
#---------------------------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_4.png', width=600, height=250)
gridExtra::grid.arrange(p_before, p_after, nrow=1)
dev.off()

