#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=200000
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice


cd 21_lassosum_imp/

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate lassoSum

Rscript lasso_script.R --use-conda

