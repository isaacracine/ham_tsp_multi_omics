#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=200000
#SBATCH -A llcg_ad_sao
#SBATCH --cluster=wice


cd /staging/leuven/stg_00092/HAM_TSP_Isaac/samples_with_exceptions/qc_feb_24_2023/21_lassosum_imp/

export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
module purge
eval "$(conda shell.bash hook)"
conda activate lassoSum

Rscript lasso_script_pvl_ref.R --use-conda

