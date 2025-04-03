#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=320G
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice
#SBATCH --ntasks-per-node=8
#SBATCH --partition=bigmem

cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/35_regenie/

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate regenie_env

regenie \
  --step 2 \
  --bed input/htlv_imp \
  --covarFile input/gwas_covs.txt  \
  --phenoFile input/pheno.txt \
  --bsize 200 \
  --bt \
  --firth --approx \
  --ref-first \
  --threads 8 \
  --strict \
  --pred output/fit_bin_out_step1_pred.list \
  --out output/test_fin_out_firth
