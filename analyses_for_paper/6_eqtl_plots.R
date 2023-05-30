rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(ggrepel)
library(MetBrewer)
library(ggpubr)
library(scales)
library(dplyr)


#---------#
# LISTS   #
#---------#

pheno_codes = c("urea_scaled", "chloride_scaled",
                "crp_scaled", "il6_scaled", "neutrophil_scaled", 
                "lymphocyte_scaled", "neutro_lympho_ratio_scaled",
                "ddimer_scaled", 'out_ICU')

pheno_txt = c('Urea', 'Chloride', 
              'C-reactive protein', 'Interleukin-6', 'Neutrophil count', 
              'Lymphocyte count', 'Neutrophil-to-lymphocyte ratio', 
              'D-dimers', 'ICU admission')

pheno_ann_df = data.frame(cbind(pheno_codes, pheno_txt))
names(pheno_ann_df) = c('Codes', 'Text')

#------------#
# FUNCTIONS  #
#------------#

merge_tables = function(table_outcomes, table_biomarkers){
  
  # Prepare tables for merging 
  table_outcomes = table_outcomes %>% select(Dependent.variable, miRNA, Effect, P.value, Significance)
  table_biomarkers = table_biomarkers %>% select(Dependent.variable, miRNA, Effect, P.value, Significance)
  
  # Merge all results in one table
  tables = rbind(table_outcomes, table_biomarkers)
  
  # Return table
  return(tables)
  
}

extract_eqtl_sumstats = function(mirna, eqtl_dir){
  
  eqtl_path = str_interp('${eqtl_dir}/${mirna}_cis_linear_covar.assoc.linear.adjusted') 
  
  # If file exists, there are some distal eQTLs. 
  if (file.exists(eqtl_path)){
    
    df = data.table::fread(eqtl_path, header=TRUE)  
    
    n_tested = nrow(df)
    
    n_sig_bonferoni = nrow(df %>% filter(BONF < 0.05))
    
    result = c(mirna, n_tested, n_sig_bonferoni)
  } else {
    result = c(mirna, 0, NA)
  }
  
  return(result)
  
}

return_combined_eqtl_df = function(mirna, eqtl_dir){
  
  eqtl_path_adjusted = str_interp('${eqtl_dir}/${mirna}_cis_linear_covar.assoc.linear.adjusted') 
  eqtl_path_unadjusted = str_interp('${eqtl_dir}/${mirna}_cis_linear_covar.assoc.linear') 
  
  if (file.exists(eqtl_path_adjusted)){
    
    # Load adjusted and unjadjusted 
    df_adjusted = data.table::fread(eqtl_path_adjusted, header=TRUE)    
    df_unadjusted = data.table::fread(eqtl_path_unadjusted, header=TRUE)    
    
    # Merge the two tables 
    df = merge(df_unadjusted %>% dplyr::select(CHR, SNP, BP, A1, NMISS, BETA, P),
               df_adjusted %>% dplyr::select(CHR, SNP, BONF, FDR_BH),
               by=c('CHR', 'SNP'))
    df$MIRNA = mirna
    
  } else {
    df = data.frame(matrix(nrow=0, ncol=10))
    names(df) = c("CHR", "SNP", "BP", "A1", "NMISS", "BETA", "P", "BONF", "FDR_BH", "MIRNA")
  }
  
  return(df)
  
}

clean_oasis_ann = function(oasis_ann){
  
  # Chromosome
  for (i in 1:nrow(oasis_ann)){
    oasis_ann[i,'Chr'] = unlist(strsplit(oasis_ann[i,'Chromo'], 'chr'))[2]
  }
  
  oasis_ann$miRNA_ID = str_replace_all(oasis_ann$miRNA_ID, '\\.', '-')
  
  return(oasis_ann)
}

annotate_main_cis_eqtl_w_phenotypes = function(tables, main_cis_eqtl_ann){
  
  covid_pheno = unique(tables$Dependent.variable)
  n_pheno = length(covid_pheno)
  
  for (pheno in covid_pheno){
    
    tmp_df = tables %>%
      filter(Dependent.variable == pheno)
    
    for (i in 1:nrow(tmp_df)){
      
      mirna = tmp_df[i,'miRNA']
      effect = tmp_df[i,'Effect']
      
      if (effect < 0){
        eff_sign = -1
      } else if (effect > 0){
        eff_sign = 1
      }
      
      main_cis_eqtl_ann[main_cis_eqtl_ann$MIRNA == mirna, pheno] = eff_sign
      
    }
  }
  
  main_cis_eqtl_ann$n_na = rowSums(is.na(main_cis_eqtl_ann %>%
                                           dplyr::select(all_of(covid_pheno))))
  
  main_cis_eqtl_ann[main_cis_eqtl_ann$n_na == n_pheno, 'associated_with_pheno'] = 'No'
  main_cis_eqtl_ann[main_cis_eqtl_ann$n_na != n_pheno, 'associated_with_pheno'] = 'Yes'
  
  return(main_cis_eqtl_ann)
  
}

annotate_pos_and_neg_associations = function(top_eqtls, pheno_ann_df){
  
  top_eqtls$Negative_associations = NA 
  top_eqtls$Positive_associations = NA 
  
  tmp_df = top_eqtls %>% 
    dplyr::select(urea_scaled:out_ICU)
  
  for (i in 1:nrow(tmp_df)){
    
    neg_assoc = c()
    pos_assoc = c()
    
    for (j in 1:ncol(tmp_df)){
      
      pheno = names(tmp_df)[j]
      eff = tmp_df[i,pheno]
      
      if (is.na(eff) == FALSE){
        if (eff < 0){
          neg_assoc = c(neg_assoc, pheno)
        } else if (eff > 0){
          pos_assoc = c(pos_assoc, pheno)
        }
      }
    }
    
    neg_assoc_txt = convert_codes_to_txt(neg_assoc, pheno_ann_df)
    pos_assoc_txt = convert_codes_to_txt(pos_assoc, pheno_ann_df)
    
    top_eqtls[i, 'Negative_associations'] = neg_assoc_txt 
    top_eqtls[i, 'Positive_associations'] = pos_assoc_txt 
    
  }
  
  return(top_eqtls)
}

convert_codes_to_txt = function(list, pheno_ann_df){
  
  len = length(list)
  
  if (len != 0){
    for (i in 1:len){
      code = list[i]
      txt = pheno_ann_df[which(pheno_ann_df$Codes == code), 'Text']
      list[i] = txt
    }
    
    list = paste(list, collapse=', ')
    
  } else {
    list = NA 
  }
  
  return(list)
}

make_snp_position = function(top_eqtls){
  
  for (i in 1:nrow(top_eqtls)){
    chr = top_eqtls[i,'CHR']
    bp = top_eqtls[i,'BP']
    
    snp_pos = str_interp("Chr${chr}:${bp}")
    top_eqtls[i,'BP'] = snp_pos
  }
  
  return(top_eqtls)
  
}

annotate_with_wannovar = function(top_eqtls, wannovar_ann, freq_df, bim){
  
  wannovar_annotated = merge(wannovar_ann, bim, 
                             by.x=c('Chr', 'Start', 'Ref', 'Alt'),
                             by.y=c('CHR', 'POS', 'A1', 'A2'), all.x=TRUE)
  
  wannovar_to_merge = wannovar_annotated %>%
    dplyr::select(SNP, Func.refGene, Gene.refGene) %>%
    rename(GeneFunction = Func.refGene,
           Gene = Gene.refGene)
  
  freq_to_merge = freq_df %>%
    dplyr::select(SNP, MAF) %>%
    mutate(MAF = round(MAF, 2))
  
  merged_df = merge(top_eqtls, freq_to_merge, by='SNP', all.x=TRUE)
  merged_2_df = merge(merged_df, wannovar_to_merge,
                      by='SNP', all.x=TRUE)
  
  return(merged_2_df)
}

make_manhattan = function(eqtl_results, P, log10_p_thresh, highlight_color, label, title){
  
  if (P == 'BONF'){
    names(eqtl_results)[8:10] = c('P_UNADJ', 'P', 'FDR_BH')
  } else if (P == 'FDR'){
    names(eqtl_results)[8:10] = c('P_UNADJ', 'P_BONF', 'P')
  }
  
  don = eqtl_results %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(eqtl_results, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(associated_with_pheno == 'Yes', "yes", "no")) %>%
    mutate( is_annotate=ifelse(-log10(P)>log10_p_thresh & associated_with_pheno == 'Yes' &
                                 top_snp == 'yes', "yes", "no")) 
  
  max_y = ceiling(max(-log10(don$P)))
  
  axisdf = don %>% 
    group_by(CHR) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  plot = ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point(color = 'grey', alpha=0.8, size=2) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y)) +     # remove space between plot area and x axis 
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size=14, hjust=0.5, vjust=0.5, face='bold')
    ) +
    
    xlab("") + ylab(expression('-log'[10]*'P')) +
    
    ggtitle(title) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "dark grey") +
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color=highlight_color, size=2)
  
  # Add label using ggrepel to avoid overlapping
  if (label == 'label_mirna'){
    plot + geom_label_repel( data=subset(don, is_annotate=="yes"), 
                             aes(label=MIRNA),
                             box.padding   = 0.35, 
                             point.padding = 0.5,
                             segment.color = 'black',
                             max.overlaps = Inf,
                             min.segment.length = 0)
    
  } else if (label == 'label_snp'){
    plot + geom_label_repel( data=subset(don, is_annotate=="yes"), 
                             aes(label=SNP),
                             box.padding   = 0.35, 
                             point.padding = 0.5,
                             segment.color = 'black',
                             max.overlaps = Inf,
                             min.segment.length = 0)
  } else {
    plot
  }
  
}


## Anotate genotyping data 
make_snp_list = function(plink){
  
  snps = data.frame(names(plink[7:ncol(plink)]))
  names(snps) = 'SNP'
  
  snps$Allele = unlist(lapply( strsplit(snps[[1]],split="_"), "[", 2))
  snps$SNP = unlist(lapply( strsplit(snps[[1]],split="_"), "[", 1))
  
  return(snps)
}

extract_barcode_from_timepoint = function(geno_ann){
  
  geno_ann$Barcode = unlist(lapply( strsplit(geno_ann[[2]],split="_"), "[", 1))
  geno_ann[which(geno_ann$Barcode == 'NGS0082-T2'), 'Barcode'] = 'NGS0082'
  
  return(geno_ann)
  
}

combine_geno_and_mirna_expression = function(plink, geno_ann, mirna_data){
  
  geno_ann = extract_barcode_from_timepoint(geno_ann)
  
  plink_ann = merge(plink, geno_ann, by.x='FID', by.y='DNA.Barcode')
  geno_mirna_df = merge(plink_ann, mirna_data, by='Barcode', all=TRUE)
  
  return(geno_mirna_df)
  
}

annotate_top_snps = function(main_cis_eqtl_ann){
  
  unique_mirna = unique(main_cis_eqtl_ann$MIRNA)
  main_cis_eqtl_ann$top_snp = 'no'
  
  for (mirna in unique_mirna){
    
    tmp_df = main_cis_eqtl_ann %>%
      filter(MIRNA == mirna) %>%
      filter(associated_with_pheno == 'Yes')
    
    if (length(unique(tmp_df$BONF)) == 1){
      if (unique(tmp_df$BONF) == 1){
        top_snp = NA
      }
    } else {
      top_snp_df = tmp_df %>%
        filter(BONF == min(tmp_df$BONF))
      
      if (nrow(top_snp_df) == 1){
        top_snp = top_snp_df$SNP
      } else {
        top_snp = top_snp_df$SNP[1]
      }
      
      i_row = which(main_cis_eqtl_ann$MIRNA == mirna &
                      main_cis_eqtl_ann$SNP == top_snp)
      
      main_cis_eqtl_ann$top_snp[i_row] = 'yes'
    }
  }
  
  return(main_cis_eqtl_ann)
}

calculate_dist_snp_to_mirna = function(main_cis_eqtl_ann){
  
  main_cis_eqtl_ann$dist_snp_start = 
    abs(as.numeric(main_cis_eqtl_ann$BP) - 
          as.numeric(main_cis_eqtl_ann$Start))
  
  main_cis_eqtl_ann$dist_snp_end = 
    abs(as.numeric(main_cis_eqtl_ann$BP) - 
          as.numeric(main_cis_eqtl_ann$End))
  
  for (i in 1:nrow(main_cis_eqtl_ann)){
    main_cis_eqtl_ann[i,'dist_snp_mirna'] =
      min(main_cis_eqtl_ann[i,'dist_snp_start'], 
          main_cis_eqtl_ann[i,'dist_snp_end'])
  }
  
  main_cis_eqtl_ann = main_cis_eqtl_ann %>%
    dplyr::select(-dist_snp_start, -dist_snp_end)
  
  return(main_cis_eqtl_ann)
  
}


## Plots for eQTLs 
plot_eqtl = function(mirna, snp, geno_mirna_df, snps_df){
  
  snp_a1 = snps_df[which(snps_df$SNP == snp), 'Allele']
  
  # Subset relevant data
  tmp_df = na.omit(geno_mirna_df[,c('Barcode', mirna, snp)])
  names(tmp_df) = c('Barcode', 'miRNA expression', 'Genotype')
  
  # Plot parameters 
  colors = met.brewer("Tiepolo", 9)[c(2,4,5)]
  
  max_y = max(tmp_df$`miRNA expression`)
  min_y = min(tmp_df$`miRNA expression`)
  comparisons_geno = list(c('0', '1'), c('1', '2'), c('0', '2'))
  
  linear_model = summary(lm(`miRNA expression` ~ Genotype, data=tmp_df))
  
  intercept = linear_model$coefficients[1,1]
  slope = linear_model$coefficients[2,1]
  pval = linear_model$coefficients[2,4]
  pval_format = formatC(pval, format = "e")
  
  mirna = str_replace_all(mirna, "\\.", "-")
  mirna = str_replace_all(mirna, "hsa-", "")
  
  # Plot
  p = ggplot(tmp_df, aes(x=factor(Genotype), y=`miRNA expression`, fill=factor(Genotype))) +
    geom_violin() + geom_jitter(width=0.2, color='black') +
    stat_summary(fun = "mean", geom = "crossbar", 
                 width = 0.5, colour = "black") +
    xlab(str_interp("Number of minor alleles")) + 
    ylab("miRNA expression") +
    geom_abline(intercept = intercept, slope = slope, 
                linetype = 'dashed', color = 'black') +
    theme_classic() +
    theme(axis.title = element_text(size=13),
          axis.text = element_text(size=13),
          strip.text = element_text(size=13),
          legend.position = 'none',
          plot.title = element_text(size=14, vjust=0.5, hjust=0.5, face='bold')) +
    scale_fill_manual(values = colors) +
    ggtitle(str_interp("${snp} and ${mirna}\n")) +
    annotate(geom="text", label = str_interp("P = ${pval_format}"),
             x = 1.2, y=max_y+0.15, size=4.5) +
    # stat_compare_means(comparisons = comparisons_geno,
    #                    label = "p.signif", method = "t.test", size=5,
    #                    label.y = c(max_y+0.2, max_y+0.5, max_y+0.8)) +
    #stat_compare_means(method='anova', label.y = max_y+0.1, size=4.5,
    #                   label.x=1.8) +
    ylim(min_y*1.1, max_y+0.2)
  
  
}

main_make_fine_mapping_plot = function(mirna){
  
  mirna_chr = as.character(unique(eqtl_list_df[which(eqtl_list_df$MIRNA == mirna), 'Chr']))
  mirna_start = as.numeric(unique(eqtl_list_df[eqtl_list_df$MIRNA == mirna, 'Start']))
  mirna_end = as.numeric(unique(eqtl_list_df[eqtl_list_df$MIRNA == mirna, 'End']))
  
  mirna_mean_pos = (mirna_start + mirna_end)/2
  
  mirna_gwas = data.frame(main_cis_eqtl_ann %>%
                            filter(MIRNA == mirna))
  
  cis_distance = 300000
  
  p = make_fine_mapping_plot(mirna, mirna_gwas, 
                             mirna_mean_pos, '#3FA796', cis_distance)
  
  return(p)
}

make_fine_mapping_plot = function(mirna, mirna_gwas, 
                                  mirna_mean_pos, highlight_color, cis_distance){
  
  don = mirna_gwas %>% 
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(associated_with_pheno == 'Yes' &
                                  BONF < 0.05, "yes", "no")) %>%
    mutate( is_annotate=ifelse(top_snp == 'yes', "yes", "no")) 
  
  max_y = ceiling(max(-log10(don$BONF)))
  range_x = max(don$BP) - min(don$BP)
  
  mirna_pos = round(mirna_mean_pos, -5)
  x_min = mirna_pos - 200000
  x_max = mirna_pos + 200000
  
  plot = ggplot(don, aes(x=BP, y=-log10(BONF))) +
    
    # Show all points
    geom_point(alpha=0.8, size=3) +
    geom_point(aes(x=mirna_mean_pos, y=-(max_y/12)), 
               color="#C74B50", size=5, shape=18) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size=14, hjust=0.5, vjust=0.5, face='bold'),
      axis.text = element_text(size=12),
      axis.title = element_text(size=13)
    ) +
    
    xlab("") + ylab(expression('-log'[10]*'P')) +
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color=highlight_color, size=2) +
    
    # Add label for top SNP 
    geom_label_repel( data=subset(don, is_annotate=="yes"), 
                      aes(label=SNP), 
                      box.padding   = 0.35, 
                      point.padding = 0.5,
                      segment.color = 'black',
                      max.overlaps = Inf,
                      min.segment.length = 0) +
    
    ggtitle(mirna) +
    scale_x_continuous(labels = comma, 
                       breaks = c(x_min, x_max))
  
  return(plot)
}


# wANNOVAR
save_wannovar_file_from_bim = function(main_cis_eqtl_ann, bim, out_filepath){
  
  mirna_to_annotate_df = main_cis_eqtl_ann %>%
    filter(BONF < 0.05) %>%
    filter(associated_with_pheno == 'Yes')
  
  mirna_to_annotate = unique(mirna_to_annotate_df$MIRNA)
  
  snps_to_annotate_df = main_cis_eqtl_ann %>%
    filter(MIRNA %in% mirna_to_annotate)
  
  snps_to_annotate = unique(snps_to_annotate_df$SNP)
  
  snps_to_annotate_bim = bim %>%
    filter(SNP %in% snps_to_annotate)
  
  wannovar_format = snps_to_annotate_bim[,c('CHR', 'POS', 'POS', 'A1', 'A2', 'SNP')]
  
  # Save file 
  write.table(wannovar_format, out_filepath, sep='\t',
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  
}



#-----------#
# LOAD DATA # 
#-----------#

# eQTL folders 
dir_eqtl_covar = "~/Desktop/miRNA/eqtl_oasis/covar/maf_5p_cis_fixed/"

# Tables with significant associations 
tables_bloodpheno = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_6.csv')
tables_icu = read.csv('~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_4A.csv')
tables = rbind(tables_bloodpheno %>% 
                 dplyr::select(Dependent.variable, miRNA, Effect, P.value), 
               tables_icu %>%
                 dplyr::select(Dependent.variable, miRNA, Effect, P.value))

# Annotations 
oasis_ann = read.delim2('~/Desktop/miRNA/oasis/miRNA_coordinates_hsa.gff3_release21_2014.txt')
geno_ann = read.csv('~/Desktop/miRNA/geno/Tracker_COVID_Genotyping.csv')
wannovar_ann = read.csv('~/Downloads/query.output.genome_summary (2).csv')

mirna_list = read.csv('~/Downloads/oasis_ann_fixed.csv', header=FALSE)

# miRNA expression data
mirna_data = read.csv('~/Desktop/miRNA/mirna_processed/n96/mirna_oasis_5_in_50_percent_standard_clinical.csv')

# MAF
freq_df = data.table::fread('~/Desktop/miRNA/eqtl/mirna_91_covid_mirna_freq.frq')

# BIM file
bim = data.table::fread('~/Downloads/covid_sample_geno01_mind01_maf05_hwe005.bim')
names(bim) = c('CHR', 'SNP', 'CM', 'POS', 'A1', 'A2')

# PLINK data
plink = data.frame(data.table::fread('~/Desktop/miRNA/geno/mirna_91_recoded_minor_allele_counts.raw'))
plink = plink %>%
  dplyr::select(!(ends_with("_HET")))


#--------------------------------------------------------#
# TABLE: Number of cis and distal eQTLs for each miRNA   # 
#--------------------------------------------------------#
summary_df = data.frame(matrix(nrow=0, ncol=3))
names(summary_df) = c('miRNA', 'N_tested', 'N_sig_Bonferoni')

# Extract cis-eQTLs and save into folders 
for (mirna in mirna_list$V1){
  summary_df[nrow(summary_df)+1, ] = c(extract_eqtl_sumstats(mirna, dir_eqtl_covar))
}

summary_df$N_tested = as.numeric(summary_df$N_tested)
summary_df$N_sig_Bonferoni = as.numeric(summary_df$N_sig_Bonferoni)

write.csv(summary_df, '~/Desktop/miRNA/eqtl_oasis/summary_cis_eqtl_oasis.csv')


n_sig_eqtls = sum(as.numeric(na.omit(summary_df$N_sig_Bonferoni)))
tmp = summary_df %>% filter(N_tested != 0)
min_tested = min(tmp$N_tested)
max_tested = max(tmp$N_tested)
n_mirna_w_no_snps = nrow(summary_df %>% filter(N_tested == 0))


print(str_interp("We tested between ${min_tested} - ${max_tested} SNPs."))

print(str_interp("There were ${n_mirna_w_no_snps} miRNAs with no tested SNPs. "))

print(str_interp("There are a total of ${n_sig_eqtls} significant cis-eQTLs in this dataset."))


#-------------------------------------------#
# DATA FRAME: All cis miRNA + annotations   #
#-------------------------------------------#
main_cis_eqtl_df = data.frame(matrix(nrow=0, ncol=10))
names(main_cis_eqtl_df) = c("CHR", "SNP", "BP", "A1", "NMISS", "BETA", "P", "BONF", "FDR_BH", "MIRNA")

# Combine all eQTLs across miRNAs 
for (mirna in mirna_list$V1){
  df_mirna = return_combined_eqtl_df(mirna, dir_eqtl_covar)
  main_cis_eqtl_df = rbind(main_cis_eqtl_df, df_mirna)
}

main_cis_eqtl_df$MIRNA = str_replace_all(main_cis_eqtl_df$MIRNA, "\\.", "-")

# Combine that with OASIS annotations
oasis_ann = clean_oasis_ann(oasis_ann)
main_cis_eqtl_ann = data.frame(merge(main_cis_eqtl_df,
                                     oasis_ann %>% 
                                       dplyr::select(miRNA_ID, Chr, Start, End), 
                                     by.x='MIRNA', by.y='miRNA_ID', allow.cartesian=TRUE))

# Anotate information about phentoype associations 
main_cis_eqtl_ann = annotate_main_cis_eqtl_w_phenotypes(tables, main_cis_eqtl_ann)
main_cis_eqtl_ann = annotate_top_snps(main_cis_eqtl_ann)
main_cis_eqtl_ann = calculate_dist_snp_to_mirna(main_cis_eqtl_ann)

tmp = main_cis_eqtl_ann %>%
  filter(BONF < 0.05) %>%
  filter(associated_with_pheno == 'Yes') %>%
  dplyr::select(SNP, MIRNA) %>%
  distinct()

print(str_interp("There are ${nrow(tmp)} significant eQTLs associated with ICU or a blood phenotype."))

tmp_within_10k_bp = main_cis_eqtl_ann %>%
  filter(BONF < 0.05) %>%
  filter(associated_with_pheno == 'Yes') %>%
  filter(dist_snp_mirna < 10000) %>%
  dplyr::select(SNP, MIRNA) %>%
  distinct()

print(str_interp("Of those, ${nrow(tmp_within_10k_bp)} were found within a 10,000 base pairs window."))

tmp_within_10k_bp = main_cis_eqtl_ann %>%
  filter(BONF < 0.05) %>%
  filter(associated_with_pheno == 'Yes') %>%
  filter(dist_snp_mirna < 10000)

View(tmp_within_10k_bp)


#-------------------------#
# TABLE 2: Top cis-eQTLs  #
#-------------------------#

# Prepare file for wANNOVAR annotations
# save_wannovar_file_from_bim(main_cis_eqtl_ann, bim,
#                             '~/Downloads/wannovar_file.txt')

top_eqtls = data.frame(main_cis_eqtl_ann %>%
                         filter(BONF < 0.05) %>%
                         filter(associated_with_pheno == 'Yes') %>%
                         filter(top_snp == 'yes') %>%
                         dplyr::select(MIRNA, CHR, SNP, BP, A1, BONF, dist_snp_mirna, urea_scaled:out_ICU))

top_eqtls = annotate_pos_and_neg_associations(top_eqtls, pheno_ann_df)
top_eqtls = make_snp_position(top_eqtls)

top_eqtls = top_eqtls %>%
  dplyr::select(MIRNA:dist_snp_mirna, Negative_associations, Positive_associations) %>%
  rename(P_BONFERONNI = BONF) %>%
  rename(SNP_TO_MIRNA_DISTANCE = dist_snp_mirna)

top_eqtl_ann = annotate_with_wannovar(top_eqtls, wannovar_ann, freq_df, bim)
top_eqtl_ann = top_eqtl_ann %>%
  arrange(P_BONFERONNI)

top_eqtl_ann[is.na(top_eqtl_ann)] = " "
top_eqtl_ann$Gene = str_replace(top_eqtl_ann$Gene, ";", ", ")

write.csv(top_eqtl_ann, '~/Desktop/covid_paper_outputs_oasis/Table_2.csv', row.names=FALSE)


#----------------------------------------------------------------#
# SUPPLEMENTARY TABLE 11: Full list of COVID-associated eQTLs    #
#----------------------------------------------------------------#
eqtls = data.frame(main_cis_eqtl_ann %>%
                     filter(BONF < 0.05) %>%
                     filter(associated_with_pheno == 'Yes') %>%
                     dplyr::select(MIRNA, SNP, CHR, BP, BONF, associated_with_pheno,
                                   dist_snp_mirna, urea_scaled:out_ICU))

eqtls = annotate_pos_and_neg_associations(eqtls, pheno_ann_df)
eqtls = make_snp_position(eqtls)

all_eqtls = eqtls %>%
  dplyr::select(MIRNA, SNP, CHR, BP, BONF, associated_with_pheno, 
                dist_snp_mirna, Negative_associations, Positive_associations) %>%
  arrange(BONF) %>%
  rename(P_BONFERONNI = BONF) %>%
  rename(SNP_TO_MIRNA_DISTANCE = dist_snp_mirna) %>%
  rename(SNP_POSITION = BP)

all_eqtls = annotate_with_wannovar(all_eqtls, wannovar_ann, freq_df, bim)
all_eqtls = all_eqtls %>%
  arrange(P_BONFERONNI)

all_eqtls[is.na(all_eqtls)] = " "
all_eqtls$Gene = str_replace(all_eqtls$Gene, ";", ", ")

write.csv(all_eqtls, '~/Desktop/covid_paper_outputs_oasis/Supplementary_Table_11.csv', row.names=FALSE)


#-----------------------------------#
# MANHATTAN PLOTS (panels A and B)  #
#-----------------------------------#
df_for_plot = main_cis_eqtl_ann
df_for_plot$MIRNA = str_replace_all(df_for_plot$MIRNA, 'hsa-', '')

p_A = make_manhattan(df_for_plot, 'BONF', -log10(0.0005), "#C74B50", 'label_mirna', '')


#--------------#
# eQTL PLOTS   #
#--------------#

## Prepare data 
snps_df = make_snp_list(plink)
names(plink)[7:ncol(plink)] = snps_df$SNP

geno_mirna_df = combine_geno_and_mirna_expression(plink, geno_ann, mirna_data)

# Select highly significant + relevant eQTLs 
eqtl_list_df = main_cis_eqtl_ann %>%
  filter(BONF < 0.05) %>%
  filter(associated_with_pheno == 'Yes') %>%
  filter(top_snp == 'yes')

eqtl_list = eqtl_list_df %>%
  dplyr::select(SNP, MIRNA) %>%
  distinct()

# Make plots for each eQTL 
for (i in 1:nrow(eqtl_list)){
  snp = as.character(eqtl_list[i, 'SNP'])
  mirna = str_replace_all(as.character(eqtl_list[i, 'MIRNA']), "-", ".")
  
  assign(str_interp("p${i}"),
         plot_eqtl(mirna, snp, geno_mirna_df, snps_df))
  
  assign(str_interp("p_${mirna}"), 
         plot_eqtl(mirna, snp, geno_mirna_df, snps_df))
}

png('~/Desktop/miRNA/eqtl_oasis/eqtls_grid.png', width=1000, height=1500)
gridExtra::grid.arrange(p1, p2, p3, p4, 
                        p5, p6, p7, p8, 
                        p9, p10, p11, p12, 
                        p13, p14, p15, 
                        p16, p17, p18, p19, p20,
                        p21, p22, p23, p24, p25,
                        p26, p27, p28, p29, p30, nrow=10)
dev.off()


#--------------------------#
# eQTL fine-mapping plots  #
#--------------------------#
count=31

for (mirna in eqtl_list$MIRNA){

  assign(str_interp("p${count}"), 
         main_make_fine_mapping_plot(mirna))
  assign(str_interp("p_${mirna}_finemap"), 
         main_make_fine_mapping_plot(mirna))
  
  count = count+1
}

png('~/Desktop/miRNA/eqtl_oasis/eqtls_fine_mapped_grid.png', width=1400, height=1600)
gridExtra::grid.arrange(p31, p32, p33, p34, p35, p36, p37, p38, p39, p40,
                        p41, p42, p43, p44, p45, p46, p47, p48, p49, p50,
                        p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, nrow=6)
dev.off()



#----------------------------------#
# eQTL violin + fine-mapping plot  #
#----------------------------------#

for (i in 1:nrow(eqtl_list)){
  snp = as.character(eqtl_list[i, 'SNP'])
  mirna = as.character(eqtl_list[i, 'MIRNA'])
  
  mirna_dot = str_replace_all(mirna, '-', '.')
  
  p_violin = plot_eqtl(mirna_dot, snp, geno_mirna_df, snps_df)
  p_fine_map = main_make_fine_mapping_plot(mirna)
  
  png(str_interp('~/Desktop/miRNA/eqtl_oasis/eqtls_${mirna}.png'), 
      width=700, height=350)
  gridExtra::grid.arrange(p_violin, p_fine_map, nrow=1)
  dev.off()
}



#-------------------------------------#
# FIGURE 3:                           #
#  (A) Manhattan plots for eQTLs      #
#  (B-C) eQTL MIR5189                 #
#  (D-E) eQTL MIR3934                 #
#  (F-G) eQTL MIR181A1                #
#-------------------------------------#

p_MIR5189_combo = ggarrange(p_hsa.miR.5189.5p, `p_hsa-miR-5189-5p_finemap`, 
                            nrow=1, labels = c('B', 'C'),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR625_combo = ggarrange(p_hsa.miR.625.3p, `p_hsa-miR-625-3p_finemap`, 
                            nrow=1, labels = c('D', 'E'),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR181A1_combo = ggarrange(p_hsa.miR.181a.3p, `p_hsa-miR-181a-3p_finemap`,
                             nrow=1, labels = c('F', 'G'),
                             label.x = 0, label.y = 1, 
                             font.label = list(size=18, color='black', face='bold'),
                             widths = c(0.4, 0.6))

p_full = ggarrange(p_A, p_MIR5189_combo, p_MIR625_combo, p_MIR181A1_combo,
                   nrow=4, labels = c('A', NA, NA, NA),
                   label.x = 0, label.y = 1, 
                   font.label = list(size=18, color='black', face='bold'))

png('~/Desktop/covid_paper_outputs_oasis/Figure_3.png', width=650, height=800)
plot(p_full)
dev.off()


#----------------------------#
# SUPPLEMENTARY FIGURE 15:   #
# eQTLs and fine-mapping     #
#----------------------------#

p_MIR7854_combo = ggarrange(p_hsa.miR.7854.3p, `p_hsa-miR-7854-3p_finemap`,
                            nrow=1, labels = c('A', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR1255_combo = ggarrange(p_hsa.miR.1255a, `p_hsa-miR-1255a_finemap`,
                            nrow=1, labels = c('B', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR4742_combo = ggarrange(p_hsa.miR.4742.3p , `p_hsa-miR-4742-3p_finemap`,
                           nrow=1, labels = c('C', NA),
                           label.x = 0, label.y = 1, 
                           font.label = list(size=18, color='black', face='bold'),
                           widths = c(0.4, 0.6))

p_MIR342_combo = ggarrange(p_hsa.miR.342.5p, `p_hsa-miR-342-5p_finemap`, 
                           nrow=1, labels = c('D', NA),
                           label.x = 0, label.y = 1, 
                           font.label = list(size=18, color='black', face='bold'),
                           widths = c(0.4, 0.6))

p_MIR339_combo = ggarrange(p_hsa.miR.339.3p, `p_hsa-miR-339-3p_finemap`, 
                           nrow=1, labels = c('E', NA),
                           label.x = 0, label.y = 1, 
                           font.label = list(size=18, color='black', face='bold'),
                           widths = c(0.4, 0.6))

p_MIR5003_combo = ggarrange(p_hsa.miR.5003.3p, `p_hsa-miR-5003-3p_finemap`, 
                            nrow=1, labels = c('F', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_7976_combo = ggarrange(p_hsa.miR.7976, `p_hsa-miR-7854-3p_finemap`, 
                         nrow=1, labels = c('G', NA),
                         label.x = 0, label.y = 1, 
                         font.label = list(size=18, color='black', face='bold'),
                         widths = c(0.4, 0.6))

p_MIR3177_combo = ggarrange(p_hsa.miR.3177.3p, `p_hsa-miR-3177-3p_finemap`,
                            nrow=1, labels = c('H', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR6777_combo = ggarrange(p_hsa.miR.6777.3p, `p_hsa-miR-6777-3p_finemap`,
                            nrow=1, labels = c('I', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR7155_combo = ggarrange(p_hsa.miR.7155.5p, `p_hsa-miR-7155-5p_finemap`,
                            nrow=1, labels = c('J', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR6868_combo = ggarrange(p_hsa.miR.6868.3p, `p_hsa-miR-6868-3p_finemap`,
                            nrow=1, labels = c('K', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_MIR4999_combo = ggarrange(p_hsa.miR.4999.5p, `p_hsa-miR-4999-5p_finemap`,
                            nrow=1, labels = c('L', NA),
                            label.x = 0, label.y = 1, 
                            font.label = list(size=18, color='black', face='bold'),
                            widths = c(0.4, 0.6))

p_total = ggarrange(p_MIR7854_combo, p_MIR1255_combo, p_MIR4742_combo, p_MIR342_combo,
                    p_MIR339_combo, p_MIR5003_combo, p_7976_combo, p_MIR3177_combo, 
                    p_MIR6777_combo, p_MIR7155_combo, p_MIR6868_combo, p_MIR4999_combo,
                    ncol=2, nrow=6)

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_15.png', width=1200, height=1500)
plot(p_total)
dev.off()


#--------------------------------#
# PREP FOR MEDIATION ANALYSIS    #
#--------------------------------#
# Subset the SNPs and miRNAs from the top_eqtls 
all_eqtls$MIRNA = str_replace_all(all_eqtls$MIRNA, "-", ".")

relevant_geno_mirna_df = geno_mirna_df %>%
  dplyr::select(Barcode, all_of(all_eqtls$SNP), all_of(all_eqtls$MIRNA))

write.csv(relevant_geno_mirna_df, '~/Desktop/miRNA/data/relevant_geno_mirna_oasis.csv', row.names=FALSE)
