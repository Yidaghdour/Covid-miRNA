#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
#SBATCH -o 2_calculate_pca_with_1000g.%J.out
#SBATCH -e 2_calculate_pca_with_1000g.%J.err
#SBATCH --mem=90G

## This script was run on the login node line by line; not sure if it works otherwise.
dir_geno_processed='/scratch/tg1407/COVID_miRNA/data/genotyping_processed/'
dir_pca_output='/scratch/tg1407/COVID_miRNA/output/pca/'
dir_1000g='/scratch/tg1407/1000Genome/data'

/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_geno_processed}/covid_sample_geno01_mind01_maf05_hwe005 \
    --bmerge ${dir_1000g}/1000Genome_commonSNPs \
    --make-bed \
    --out ${dir_geno_processed}/mirna_and_1000g

# Identify SNPs to remove (multiple positions and 3+ alleles )
/scratch/tg1407/R/bin/Rscript extract_multiple_pos_snps.R

sort ${dir_geno_processed}/covid_sample_geno01_mind01_maf05_hwe005.bim \
     ${dir_1000g}/1000Genome_commonSNPs.bim|uniq -d > common_lines.txt
awk '{print $2}' common_lines.txt > common_snps.txt

# Remove SNPs with 3+ alleles from both datasets
/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_geno_processed}/covid_sample_geno01_mind01_maf05_hwe005 \
    --extract common_snps.txt \
    --exclude ${dir_geno_processed}/mirna_and_1000g_snps_to_remove.txt \
    --make-bed \
    --out ${dir_geno_processed}/mirna_91_removed_snps

/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_1000g}/1000Genome_commonSNPs \
    --extract common_snps.txt \
    --exclude ${dir_geno_processed}/mirna_and_1000g_snps_to_remove.txt \
    --make-bed \
    --out ${dir_1000g}/1000Genome_commonSNPs_removed_snps

# Re-merge
/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_geno_processed}/mirna_91_removed_snps \
    --bmerge ${dir_1000g}/1000Genome_commonSNPs_removed_snps \
    --make-bed \
    --out ${dir_geno_processed}/mirna_and_1000g

/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_geno_processed}/mirna_and_1000g \
    --pca \
    --out ${dir_pca_output}/mirna_and_1000g_pca

rm common_snps.txt
rm common_lines.txt
rm mirna_*
