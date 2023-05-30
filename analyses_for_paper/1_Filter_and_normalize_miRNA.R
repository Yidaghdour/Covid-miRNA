rm(list=ls())


#--------------#
# LIBRARIES    #
#--------------#
library(dplyr)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(ggpubr)

library(devtools)
library(factoextra)
library(stringr)


#--------------#
# FUNCTIONS    #
#--------------#
make_dir_for_each_method = function(methods){
  for (method in methods){
    
    # Make directory for method 
    method_dir_name = paste("~/Desktop/miRNA/", method, sep='')
    dir.create(method_dir_name, showWarnings = FALSE)
    
    # Make directory for qc_plots 
    dir_name = paste(method_dir_name, "/qc_plots", sep='')
    dir.create(dir_name, showWarnings=FALSE)
  }
}

set_wd_to_qc_plots = function(method){
  dir_name = paste("~/Desktop/miRNA/", method, "/qc_plots", sep='')
  setwd(dir_name)
}


# Supplementary functions for load_and_process_data()
read_mirna_data = function(method){
  if (method == "encode"){
    df = readxl::read_xlsx("~/Desktop/miRNA/encode/Consolidated_Count.txt_Annotated.xlsx")
  } else if (method == "cgsb"){
    df = read.csv("~/Desktop/miRNA/cgsb/miRNA_output.csv")
  } else if (method == 'oasis'){
    df = read.delim2("~/Desktop/miRNA/oasis/Consolidated_Count.txt")
  }
  
  return(df)
}

format_encode = function(df){
  
  # Get names of miRNa's 
  mirna_names = df$GeneName
  
  # Remove unnecessary columns 
  cols_to_drop = c("GeneName", "Location", "GeneID")
  df = df[, !(names(df) %in% cols_to_drop)]
  
  # Transpose data: colnames = individual IDs, rownames = miRNA's 
  df_t = data.frame(t(df))
  colnames(df_t) = mirna_names 
  
  # Print outputs
  n_individuals = nrow(df_t)
  n_mirnas = ncol(df_t)
  
  print(paste("There are a total of", n_individuals, "individuals."))
  print(paste("and", n_mirnas, "unqiue miRNA's before filtering."))
  
  return(df_t)
}

format_cgsb = function(df){
  
  # Get names of miRNa's 
  mirna_names = df$GeneID
  
  # Remove unnecessary columns 
  cols_to_drop = c("GeneID")
  df = df[, !(names(df) %in% cols_to_drop)]
  
  # Transpose data: colnames = individual IDs, rownames = miRNA's 
  df_t = data.frame(t(df))
  colnames(df_t) = mirna_names 
  
  # Print outputs
  n_individuals = nrow(df_t)
  n_mirnas = ncol(df_t)
  
  print(paste("There are a total of", n_individuals, "individuals."))
  print(paste("and", n_mirnas, "unqiue miRNA's before filtering."))
  
  return(df_t)
}

format_oasis = function(df){
  
  # Get names of miRNa's 
  mirna_names = df$miRNA
  
  # Remove unnecessary columns 
  cols_to_drop = c("miRNA")
  df = df[, !(names(df) %in% cols_to_drop)]
  
  # Transpose data: colnames = individual IDs, rownames = miRNA's 
  df_t = data.frame(t(df))
  colnames(df_t) = mirna_names 
  
  # Print outputs
  n_individuals = nrow(df_t)
  n_mirnas = ncol(df_t)
  
  print(paste("There are a total of", n_individuals, "individuals."))
  print(paste("and", n_mirnas, "unqiue miRNA's before filtering."))
  
  return(df_t)
}


convert_ids_to_clinical_encode = function(df, clinical_barcodes){
  
  # Convert individual IDs to standard IDs 
  for (i in 1:nrow(df)){
    
    # Extract miRNA barcode
    barcode = rownames(df)[i]
    barcode_core = unlist(strsplit(barcode, "_"))[1]
    
    if (nchar(barcode_core) == 2){
      barcode_core = paste("NGS00", barcode_core, sep='')
    }
    
    # Find applicable clinical barcode
    for (id in clinical_barcodes){
      if (grepl(barcode_core, id, fixed=TRUE)){
        clinical_barcode = id
      }
    }
    
    # Replace 
    rownames(df)[i] = clinical_barcode
  }
  
  return(df)
  
}

convert_ids_to_clinical = function(df, clinical_barcodes){
  
  # Convert individual IDs to standard IDs 
  for (i in 1:nrow(df)){

    # Extract miRNA barcode
    barcode = rownames(df)[i]
    barcode_core = unlist(strsplit(barcode, "_"))[1]
    
    if (method %in% c('oasis', 'cgsb')){
      barcode_core = substr(barcode_core, 2, 4)
    }

    if (nchar(barcode_core) == 2){
      barcode_core = paste("NGS00", barcode_core, sep='')
    }
    
    # Find applicable clinical barcode
    for (id in clinical_barcodes){
      if (grepl(barcode_core, id, fixed=TRUE)){
        clinical_barcode = id
      }
    }
    
    # Replace 
    rownames(df)[i] = clinical_barcode
  }
  
  return(df)
}


load_and_process_data = function(method, clinical_barcodes){
  
  print(paste("Processing miRNA data obtained from", method))
  
  # Read data 
  data = read_mirna_data(method)
  
  # Clean the formatting and annotate with clinical IDs 
  if (method == 'encode'){
    data_formatted = format_encode(data)
  } else if (method == 'cgsb'){
    data_formatted = format_cgsb(data)
  } else if (method == 'oasis'){
    data_formatted = format_oasis(data)
  }
  
  data_processed = convert_ids_to_clinical(data_formatted, clinical_barcodes)
  
  return(data_processed)
}

filter_data = function(data, filter){
  n_individuals = nrow(data)
  n_mirnas = ncol(data)
  
  # Get key variables from filter 
  if (filter == "5_in_50_percent"){
    min_count = 5 
    min_indiv = n_individuals/2
  } else if (filter == '1_in_all_samples'){
    min_count = 1
    min_indiv = n_individuals
  } else if (filter == '1_in_any_sample'){
    min_count = 1 
    min_indiv = 1
  } else {
    print("There is an error. The filter doesn't match one of the three options:")
    print("minimum 5 in 50% of samples, minimum 1 in all samples, OR minimum 1 in any 1 sample")
  }
  
  # Identify miRNA's to drop, based on filter criteria 
  mirna_to_drop = c()
  
  for (j in 1:ncol(data)){
    n_indiv_over_min_count = sum(data[1:n_individuals,j] > min_count)
    
    if (n_indiv_over_min_count < min_indiv){
      mirna_name = names(data)[j]
      mirna_to_drop = c(mirna_to_drop, mirna_name)
    }
  }
  
  # Drop miRNA
  data_filtered = data[, !(names(data) %in% mirna_to_drop)]
  
  # Print summary 
  n_mirna_dropped = length(mirna_to_drop)
  n_mirna_remain = ncol(data_filtered)
  
  print(paste(n_mirna_dropped, "miRNAs were dropped based on filter", filter))
  print(paste(n_mirna_remain, "miRNAs remain after filtering."))
  
  return(data_filtered)
}

log2_transform_data = function(data){
  
  # Add +1 to the data before log-transforming
  data = data + 1
  
  # Log2 transform data 
  data = log2(data)
  
  # Transpose data 
  data = data.frame(t(data))
  
  print("Added +1 to all values, log2 transformed and transposed the matrix.")
  print("Individuals are now columns, miRNAs are rows.")
  
  return(data)
}

normalize_data = function(data, normalization, outliers, perc_missingness){

  if (normalization == 'Standard'){

    for (c in 1:ncol(data)){

      # Demean and divide by SD
      mean = mean(na.omit(data[,c]))
      sd = sd(na.omit(data[,c]))
      data[,c] = (data[,c] - mean)/sd
      
    }
    
  } else if (normalization == 'Mean'){

    # Demean 
    for (c in 1:ncol(data)){
      mean = mean(na.omit(data[,c]))
      data[,c] = data[,c] - mean 
    }
    
  } else if (normalization == 'IQR'){
    
    # Subtract the median and divide by IQR 
    for (c in 1:ncol(data)){
      median = median(na.omit(data[,c]))
      iqr = IQR(na.omit(data[,c]))
      
      data[,c] = (data[,c] - median)/iqr
    }
    
  } else if (normalization == 'Median'){

    # Subtract the median and normalize by the median deviation 
    for (c in 1:ncol(data)){
      median = median(na.omit(data[,c]))
      tmp = abs(data[,c] - median)
      median_deviation = median(na.omit(tmp))
      
      data[,c] = (data[,c] - median)/median_deviation
    }
  }
  
  if (outliers == 'remove_outliers'){
    data[is.na(data)] = NA
    data[data == 'Inf'] = NA
    
    col_name_drop = c()
    
    for (c in 1:ncol(data)){
      n_na = sum(is.na(data[,c]))
      
      perc = n_na/nrow(data)
      
      if (perc > perc_missingness){
        col_name = names(data)[c]
        col_name_drop = c(col_name_drop, col_name)
      }
    }
    
    data = data[, !(names(data) %in% col_name_drop)]
    n_dropped = length(col_name_drop)
    
    print(str_interp("${n_dropped} individuals were removed from the analysis due to large amount of missing data."))
  }
  
  print(paste("Data was normalized using", normalization, "normalization."))
  
  return(data)
  
}

filter_data_individuals = function(data, perc_missingness){
  
  res = data.frame(colSums(is.na(data)))
  names(res) = 'n_missing'
  res$perc = round((res$n_missing / nrow(data)), 2)*100
  
  ids_to_remove = rownames(res[res$perc > perc_missingness,])
  n_to_remove = length(ids_to_remove)
  
  print(str_interp("${n_to_remove} individuals were removed due to missing data for > ${perc_missingness}% miRNA."))
  
  data = data[,!(names(data) %in% ids_to_remove)]
  
  return(data)
  
}

remove_outliers = function(df, n_sd){
  
  for (i in 1:nrow(df)){
    
    data_row = as.numeric(df[i,])
    n_values_before = sum(is.na(data_row) == FALSE)
    
    mean_mirna = mean(na.omit(data_row))
    sd_mirna = sd(na.omit(data_row))
    
    for (j in 1:ncol(df)){
      value = df[i,j]
      
      if (is.na(value) == FALSE){
        if (value > mean_mirna+n_sd*sd_mirna){
          df[i,j] = NA 
        } else if (value < mean_mirna-n_sd*sd_mirna){
          df[i,j] = NA 
        }
      }
    }
    
    data_row = as.numeric(df[i,])
    n_values_after = sum(is.na(data_row) == FALSE)
    
    n_dropped = n_values_before - n_values_after
    
    print(str_interp('${n_dropped} miRNAs dropped as outliers (outside ${n_sd} SD of mean).'))
  }
  
  return(df)
}

return_and_save_pca_plot = function(data, x, y, method, filter, norm_type, color=NULL, extra=NULL){
  
  data[is.na(data)] = 0
  n = ncol(data)
  
  # Turn data frame into matrix form 
  data_matrix = as.matrix(t(data))
  
  # Compute PCA
  pca_result = prcomp(data_matrix, center = TRUE)
  
  # Title
  if (is.null(extra)){
    plot_title = paste(method, "\n", filter, "\n", norm_type, "normalization")
  } else {
    plot_title = paste(method, "\n", filter, "\n", norm_type, "normalization,", extra)
  }
    
  # Plot PCx vs PCy
  p = autoplot(pca_result, x=x, y=y) + theme_light() +
    ggtitle(plot_title) +
    theme(axis.title = element_text(size=13),
          axis.text = element_text(size=12),
          plot.title = element_text(size=14, hjust=0.5, face="bold"))
   
  if (is.null(color) == FALSE){
    p = autoplot(pca_result, x=x, y=y, color=color) + theme_light() +
      theme(axis.title = element_text(size=13),
            axis.text = element_text(size=12))
  }

  if (is.null(extra)){
    filename = paste("PC", x, "_vs_PC", y, "_", filter, 
                     "_", norm_type, "_n_", n, ".png", sep='')
  } else {
    filename = paste("PC", x, "_vs_PC", y, "_", filter, 
                     "_", norm_type, "_", extra, "_n_", n, ".png", sep='')
  }
  
  # Save plot as PNG
  png(filename, width=500, height=400)
  print(p)
  dev.off()
  
  return(p)
  
}

return_and_save_density_plot = function(data, method, filter, before_after, palette, norm_type=NULL, extra=NULL){
  
  indiv_ids = names(data)
  n_indiv = length(indiv_ids)
  
  data_long = data %>% 
    gather(indiv_ID, mirna_Expression, all_of(indiv_ids))
  
  n = length(indiv_ids)
  
  if (before_after == 'before'){
    plot_title = paste(paste("Before normalization"))
    file_name = paste("density_", filter, "_before_normalization", "_n_", n, ".png", sep='')
  } else {
    plot_title = paste(paste(norm_type, "normalization"))
    file_name = paste("density_", filter, "_", norm_type, "_normalization_", before_after, "_n_", n, 
                      ".png", sep='')
  }
  
  p = ggplot(data_long, aes(x=mirna_Expression, color=indiv_ID)) +
    geom_density() + theme(legend.position='none') + theme_light() +
    ggtitle(plot_title) + 
    theme(axis.title = element_text(size=13),
          axis.text = element_text(size=12),
          plot.title = element_text(size=14, hjust=0.5, face="bold"),
          legend.position='none') +
    xlab("\nmiRNA expression") + ylab("Density\n") + 
    scale_color_manual(values=get_palette(palette = palette, n_indiv))
  
  png(file_name, width=450, height=320)
  plot(p)
  dev.off()
  
  return(p)
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


#--------#
# LISTS  #
#--------#

methods = c('encode', 'cgsb', 'oasis')
filters = c("5_in_50_percent", "1_in_any_sample", "1_in_all_samples")
norm_types = c('Mean', 'Median', 'IQR', 'Standard')

filters_txt = c('minimum miRNA count of 5 in 50% of samples',
                'minimum miRNA count of 1 in any sample',
                'minimum miRNA count of 1 in all samples')

filters_df = data.frame(cbind(filters, filters_txt))


#--------------------------------------#
# Load clinical and viral load data    #
#--------------------------------------#

# Load clinical data 
df_clinical_viral = read.csv("~/Desktop/COVID_miRNA/Processed_data/data_merged_w_viral.csv")
df_clinical_viral = df_clinical_viral[df_clinical_viral$Remove == 'No', ]

clinical_barcodes = df_clinical_viral$Barcode


#---------------#
# PARAMETERS    #
#---------------#

# Make directories for each method, and for qc_plots in each methods folder
make_dir_for_each_method(methods)

## CLEANING ----------------------------------------------------------------
k=1

summary_df = data.frame(matrix(
  nrow=length(methods)*length(filters)*length(norm_types)*3,
  ncol=6))

method_count = 1

for (method in methods){
  
  print(str_interp("Preparing plots for ${method}"))
  
  for (filter in filters){
    
    count = 1

    for (norm_type in norm_types){
      
      outlier_filters = c(100)
      
      for (outlier_filter in outlier_filters){
        
        # 0. Set working directory to method/qc_plots
        set_wd_to_qc_plots(method)
        
        # 1. Load, clean and annotate data -------------------------------------
        data_processed = load_and_process_data(method, clinical_barcodes)
        
        # 2. Perform filtering on data -----------------------------------------
        data_filtered = filter_data(data_processed, filter)
        
        # 3. Log transform data ------------------------------------------------
        data_transformed = log2_transform_data(data_filtered)
        
        # 4. Replace 0's with NA -----------------------------------------------
        data_transformed[data_transformed == 0] = NA
        
        # 5. Remove individuals with missingness above threshold ---------------
        data_transformed_filtered = filter_data_individuals(data_transformed, outlier_filter)
        
        # 6. Normalize data ----------------------------------------------------
        data_normalized = normalize_data(data_transformed_filtered, norm_type, 'dont_remove')
        
        # 6.5 Remove outliers (outside 3 SD) -----------------------------------
        if (filter != "1_in_any_sample"){
          data_normalized = remove_outliers(data_normalized, 3)
        }
        
        # 7. Save density plot (before/after normalization) --------------------
        p_before = return_and_save_density_plot(data_transformed, method, filter, 'before', 'Oranges')
        p_after = return_and_save_density_plot(data_normalized, method, filter, 'after', 'Purples', norm_type)
        
        # 8. Save PCA ----------------------------------------------------------
        return_and_save_pca_plot(data_transformed, 1, 2, method, filter, norm_type, extra="before_normalization")
        return_and_save_pca_plot(data_normalized, 1, 2, method, filter, norm_type, extra="after_normalization")
        
        # 9. Merge with clinical data ------------------------------------------
        data_mirna = data.frame(t(data_normalized))
        data_mirna$Barcode = rownames(data_mirna)
        
        data_main = merge(data_mirna, df_clinical_viral, by='Barcode', 
                          all.x=TRUE, all.y=FALSE)
        
        # 10. Save data --------------------------------------------------------
        if (outlier_filter == 100){
          setwd('~/Desktop/miRNA/mirna_processed/n96/')
        } else if (outlier_filter == 20){
          setwd('~/Desktop/miRNA/mirna_processed/n84/')
        } else (
          setwd('~/Desktop/miRNA/mirna_processed/n79/')
        )
        
        filename = str_interp("mirna_${method}_${filter}_${norm_type}_clinical.csv")
        write.csv(data_main, filename)
        
        # 11. Add information to a summary table -------------------------------
        n_mirna = nrow(data_normalized)
        n_indiv = ncol(data_normalized)
        
        summary_df[k, ] = c(method, filter, norm_type, outlier_filter, n_mirna, n_indiv)
        k=k+1
        
        # 12. Add plots to list -------------------------------------------------
        if (count == 1){
          assign(str_interp('p${method_count}_${count}'), p_before)
          count = count + 1
          
          assign(str_interp('p${method_count}_${count}'), p_after)
          count = count + 1
          
        } else {
          assign(str_interp('p${method_count}_${count}'), p_after)
          count = count + 1
        }
        
      }
    }
  }
  
  method_count = method_count + 1
  
}

#------------------------------------------------------------#
# SUPPLEMENTARY FIGURE 2:                                    #
# Before and after normalization for 4 different methods     #
#------------------------------------------------------------#

png('~/Desktop/covid_paper_outputs_oasis/Supplementary_Figure_2.png', width=600, height=250)
gridExtra::grid.arrange(p3_1, p3_5, nrow=1)
dev.off()


