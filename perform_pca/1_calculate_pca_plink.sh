#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
#SBATCH -o 1_calculate_pca_plink.%J.out
#SBATCH -e 1_calculate_pca_plink.%J.err
#SBATCH --mem=90G

dir_geno_processed='/scratch/tg1407/COVID_miRNA/data/genotyping_processed/'
dir_pca_output='/scratch/tg1407/COVID_miRNA/output/pca/'

/scratch/mv83/Software/plink_beta5.3/plink \
    --bfile ${dir_geno_processed}/covid_sample_geno01_mind01_maf05_hwe005 \
    --no-pheno \
    --allow-no-sex \
    --pca \
    --out ${dir_pca_output}/mirna_pca
