#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice
#SBATCH --ntasks-per-node=8

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/35_regenie/

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate regenie_env

regenie \
  --step 1 \
  --bed input/htlv_unimp \
  --covarFile input/gwas_covs.txt  \
  --phenoFile input/pheno.txt \
  --bsize 100 \
  --qt \
  --ref-first \
  --extract input/unimp_alleles_with_more_than_5_counts_pvl.snplist \
  --phenoCol log_pvl \
  --threads 8 \
  --strict \
  --out output/fit_bin_out_step1_pvl
