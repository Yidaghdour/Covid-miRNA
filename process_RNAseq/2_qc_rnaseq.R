library(tidyr)
library(stringr)
library(ggplot2)
library(tidyselect)
library(dplyr)


#------------#
# FUNCTIONS  #
#------------#

log10_transform_data = function(data){

  # Add +1 to the data before log-transforming
  data = data + 1

  # Log2 transform data
  data = log10(data)

  # Transpose data
  data = data.frame(t(data))

  print("Added +1 to all values, log10 transformed and transposed the matrix.")
  print("Individuals are now columns, miRNAs are rows.")

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

remove_outliers = function(df, n_sd){

  for (i in 1:nrow(df)){

    data_row = as.numeric(df[i,])

    mean_mirna = mean(na.omit(data_row))
    sd_mirna = sd(na.omit(data_row))

    upper_bound = mean_mirna+n_sd*sd_mirna
    lower_bound = mean_mirna-n_sd*sd_mirna

    under_lower = is.na(df[i,]) == FALSE & df[i,] < lower_bound
    above_upper = is.na(df[i,]) == FALSE & df[i,] > upper_bound

    if (sum(under_lower) != 0){
      df[i,][under_lower] = NA
    }

    if (sum(above_upper) != 0){
      df[i,][above_upper] = NA
    }

  }

  return(df)
}


#-------------#
# LOAD DATA   #
#-------------#

# RNA-seq data
data_rnaseq = read.csv('/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_filtered_transcript.csv')
rownames(data_rnaseq) = trimws(data_rnaseq$transcript_id)

# Remove NGS007 and NGS0409

data_rnaseq = data_rnaseq %>%
  select(NGS0709:NGS0427) %>%
  select(-NGS007, -NGS0409)

data_rnaseq = data.frame(t(data_rnaseq))


#----------------------#
# QC OF RNA-SEQ DATA   #
#----------------------#

# 1. Log10 transform data
data_transformed = log10_transform_data(data_rnaseq)

# 2. Replace 0's with NA
data_transformed[data_transformed == 0] = NA

# 3. Filter individuals
data_transformed_filtered = filter_data_individuals(data_transformed, 50)

# 4. Normalize RNA-seq data
data_normalized = normalize_data(data_transformed_filtered, 'Standard', 'dont_remove')

# 5. Remove outliers (+/- 3 SD)
data_normalized = remove_outliers(data_normalized, 3)

# 6. Keep fewer decimals
data_normalized_smaller = round(data_normalized, 5)
data_normalized_smallest = round(data_normalized_smaller, 2)

data_transformed_smaller = round(data_transformed, 5)
data_transformed_smallest = round(data_transformed, 2)

# 7a. Save normalized data
write.csv(data_normalized,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_normalized.csv')

write.csv(data_normalized_smaller,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_normalized_smaller.csv')

write.csv(data_normalized_smallest,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_normalized_smallest.csv')

# 7b. Save transformed data
write.csv(data_transformed,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_transformed.csv')

write.csv(data_transformed_smaller,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_transformed_smaller.csv')

write.csv(data_transformed_smallest,
          '/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_data_transformed_smallest.csv')
