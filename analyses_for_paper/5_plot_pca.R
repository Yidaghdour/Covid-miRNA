rm(list=ls())

library(MetBrewer)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

#------------#
# LOAD DATA  #
#------------#
setwd('~/Desktop/miRNA/pca/')

pca = read.table('mirna_and_1000g_pca.eigenvec', sep=' ', header=FALSE)
names(pca) = c('FID', 'IID', paste('PC', 1:20, sep=''))

clinical_df = read.csv('~/Desktop/miRNA/mirna_processed/n96/mirna_encode_5_in_50_percent_standard_clinical.csv')
clinical_df = clinical_df %>%
  select(eid)

#------------------------#
# PREPARE DATA FOR PCA   #
#------------------------#
n_cov_dataset = 91

pca$Dataset = c(rep('COVID dataset', n_cov_dataset), 
                rep('1000 Genome dataset', nrow(pca) - n_cov_dataset))

pca$Population = 'COVID-19 dataset'

for (pop in c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')){
  
  ids = read.delim2(str_interp("${pop}_IDs.txt"), header=FALSE)
  names(ids) = c('FID', 'IID')
  
  pca[pca$FID %in% ids$FID, 'Population'] = pop
}

pca = rbind(pca %>% filter(Dataset == '1000 Genome dataset'),
            pca %>% filter(Dataset == 'COVID dataset'))

#------------#
# PLOT PCA   #
#------------#
colors = met.brewer("Hiroshige", 5)

pca$Population = factor(pca$Population,
                        levels = c('AFR', 'EAS', 'EUR', 'SAS', 'COVID-19 dataset'))

p1 = ggplot(pca %>%
         filter(Population != 'AMR'), 
                aes(x=PC1, y=PC2, color=Population)) +
  geom_point(size=3,) +
  theme_classic() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        plot.title = element_text(size=15, vjust=0.5, hjust=0.5, face='bold'),
        legend.position = 'none') +
  xlab('\nPC1') + ylab('PC2\n') +
  scale_color_manual(values = c("#cee397", "#E9DAC1", "#f5b461", "#ec524b", 'black'))

p2 = ggplot(pca %>%
              filter(Population != 'AMR'), 
            aes(x=PC3, y=PC4, color=Population)) +
  geom_point(size=3,) +
  theme_classic() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        plot.title = element_text(size=15, vjust=0.5, hjust=0.5, face='bold'),
        legend.position = 'right') +
  xlab('\nPC3') + ylab('PC4\n') +
  scale_color_manual(values = c("#cee397", "#E9DAC1", "#f5b461", "#ec524b", 'black'))

#------------------------------#
# SUPPLEMENTARY FIGURE 14:     #
# PCA of genotyping data       #
#------------------------------#
png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_14.png', width=800, height=300)
ggarrange(p1, p2, nrow=1, widths=c(0.4, 0.6),
          labels = c('A', 'B'),
          label.x = 0, label.y = 1,
          font.label = list(size=18, color='black', face='bold'))
dev.off()

