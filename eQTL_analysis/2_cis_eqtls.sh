#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00
#SBATCH -o cis_eqtls_fixed_multiple_locations.%J.out
#SBATCH -e cis_eqtls_fixed_multiple_locations.%J.err
#SBATCH --mem=90G

dir_geno_processed='/scratch/tg1407/COVID_miRNA/data/genotyping_processed/'
dir_eqtl_covar_output='/scratch/tg1407/COVID_miRNA/output/oasis/eqtl_covars/'
dir_mirna_pheno='/scratch/tg1407/COVID_miRNA/data/mirna_oasis/'

cat /scratch/tg1407/COVID_miRNA/data/mirna_oasis/oasis_ann_fixed.csv | while read mirna_chr_pos
do

  mirna="$(echo ${mirna_chr_pos} | cut -d',' -f1 | xargs)"
  chr=$(echo ${mirna_chr_pos} | cut -d',' -f2 | xargs)
  start_cis_pos=$(echo ${mirna_chr_pos} | cut -d',' -f3 | xargs)
  end_cis_pos=$(echo ${mirna_chr_pos} | cut -d',' -f4 | xargs)

  for i in 1 2 3 4 5
  do

    # If it is a specific isoform (e.g. has _1)
    if [[ ${mirna} == *_${i} ]]
      then
        mirna_core=${mirna::-2}
      else
        mirna_core=${mirna}
    fi

    ## eQTL analysis
    /scratch/mv83/Software/plink_beta5.3/plink \
        --bfile ${dir_geno_processed}/covid_sample_geno01_mind01_maf05_hwe005 \
        --no-pheno \
        --allow-no-sex \
        --chr ${chr} \
        --from-kb ${start_cis_pos} \
        --to-kb ${end_cis_pos} \
        --pheno ${dir_mirna_pheno}/pheno_${mirna_core}.txt \
        --covar ${dir_mirna_pheno}/covars_age_gender.txt \
        --covar-name AGE,GENDER \
        --linear hide-covar \
        --adjust \
        --out ${dir_eqtl_covar_output}/maf_5p_cis_fixed/${mirna}_cis_linear_covar

  done
done
